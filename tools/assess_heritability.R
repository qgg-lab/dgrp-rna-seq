#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 5/3/16
#
# assess_heritability.R
#
# IMPT NOTE: This script is no longer being developed.  
# Improvements to genetic variance modeling are now implemented in gen_var_model.R
#
# Goal: Load in normalized data and apply linear mixed model to assess H^2 for each gene/feature
# 
# Usage:
# Rscript assess_heritability.R [OPTIONS] EXPR=subdir/data_[rpm|fpkm|rle].txt
#  EXPR=  The file containing expression values to assess heritability on
#  SAMPLES= The table with sample information (defaults to sample_master_table.txt)
#  WOLADJ=  TRUE/FALSE (Default=TRUE), adjust for Wolbachia effects before computing H^2 and line means
#  POOLED=  TRUE/FALSE (Default=FALSE), Run a pooled model instead of within-sex models
#  TRMT=  Name of an additional column in SAMPLES to use as Treatment variable in model (default=None)
#  CORES= The number of CPUs to use for parallelization (defaults to 1)
#   NOTE: If submitting your job through SLURM, make sure to set -c option as well
#  NQ=    TRUE/FALSE (Default=FALSE), apply Normal Quantile transformation before running H^2
#  VRFILE=  File with gene x line variant rates. If blank, no adjustment for this is performed (Default)
#           But if a file is specified, then modeling will be performed ONLY on the intersect set of genes
#           Output will be given a "_VR_" designation code
#  PERM=  If specified, permute the sample labels within each group (sex)
#         This should be provided as an integer value, which is used to both seed the RNG
#         And give the output a unique ID (output will have _PermN_ designation code)
#         All output goes in Perm/ subdir of Out Dir as well
#
# Automatically puts all output in the same subdir as the data file, with similar basename
# Normally takes log2 of the data (with pseudo count), UNLESS data file ends in *_rle
# Automatically uses sample_master_table.txt for sample information
#


# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")
# setwd("~/Projects/Tanya_EtOH_Testing/")

options(stringsAsFactors=F)

# Main package for the linear mixed model
# If not installed, run:
# install.packages(c("lme4","lmerTest","doMC"))
# NOTE: Only need to load lmerTest - it loads what it needs from lme4 or has the same functionality
# library(lme4)
library(lmerTest)
# This one is needed for parallelization later on
library(doMC)


# -- Fixed Parameters -- #
# These could be made adjustable in the command-line args
pseudo.count <- 0.001
fdr.cutoff <- 0.05
line.file <- "~ljeveret/Resources/DGRP/Freeze2/line_wolbachia_inversion_info.txt"
# line.file <- "~/Resources/DGRP/Freeze2/line_wolbachia_inversion_info.txt"


# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript assess_heritability.R [OPTIONS] EXPR=subdir/data_[rpm|fpkm|rle].txt"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("EXPR=transposons/combined_transposon_filtered_rpm.txt")
# my.args <- c("EXPR=transposons/combined_transposon_filtered_rpm.txt", "POOLED=T")
# my.args <- c("EXPR=known_rep_novel_genes/combined_samples_known_novel_fpkm.txt")
# my.args <- c("EXPR=known_rep_novel_genes/combined_samples_known_novel_fpkm.txt","VRFILE=known_rep_novel_genes/known_gene_variant_rates.txt", "CORES=4")
# my.args <- c("EXPR=known_rep_novel_genes/combined_samples_known_novel_fpkm.txt","WOLADJ=F","PERM=10")
# my.args <- c("EXPR=known_rep_novel_genes/combined_samples_known_novel_fpkm.txt", "POOLED=T", "CORES=4")
# my.args <- c("EXPR=combined.samples_DGRP_E1REG.test.txt", "SAMPLES=combined.master_table_DGRP_REGE1.txt", "CORES=4", "TRMT=ETOH")
# my.args <- c("EXPR=combined.samples_DGRP_E1REG.test.txt", "SAMPLES=combined.master_table_DGRP_REGE1.txt", "CORES=4", "TRMT=ETOH", "POOLED=T")
if(length(my.args) == 0) {
  stop("Requires at least one argument. ", usageStr)
} else {
  cat("Running with params:",my.args,"\n")
}
# Parse param names and values
argSplit=strsplit(my.args, split="=", fixed=T)
my.args <- lapply(argSplit, "[", 2)
names(my.args) <- unlist(lapply(argSplit, "[", 1))

# Fill in default value for missing params
setDefault <- function(param, default) {
  if(!(param %in% names(my.args))) {
    cat("Setting ",param,"=",default," by default\n", sep="")
    my.args[[param]] <<- default
  }
}

# Make sure EXPR param exists
if(!("EXPR" %in% names(my.args))) {
  stop("Missing EXPR parameter. ", usageStr)
}

setDefault("SAMPLES","sample_master_table.txt")

setDefault("NQ",F)
# If NQ is "0" or "1" map these to "FALSE" and "TRUE"
if(is.character(my.args$NQ) & (my.args$NQ == "0")){my.args$NQ <- F}
if(is.character(my.args$NQ) & (my.args$NQ == "1")){my.args$NQ <- T}
# Convert any other type to logical
if(!is.logical(my.args$NQ)){my.args$NQ <- as.logical(my.args$NQ)}
# If conversion produced NA, just set to F
if(is.na(my.args$NQ)) {
  cat("WARNING: Did not recognize value of NQ parameter, setting to F by default.\n")
  my.args$NQ <- F
}

setDefault("WOLADJ",T)
# If WOLADJ is "0" or "1" map these to "FALSE" and "TRUE"
if(is.character(my.args$WOLADJ) & (my.args$WOLADJ == "0")){my.args$WOLADJ <- F}
if(is.character(my.args$WOLADJ) & (my.args$WOLADJ == "1")){my.args$WOLADJ <- T}
# Convert any other type to logical
if(!is.logical(my.args$WOLADJ)){my.args$WOLADJ <- as.logical(my.args$WOLADJ)}
# If conversion produced NA, just set to T
if(is.na(my.args$WOLADJ)) {
  cat("WARNING: Did not recognize value of WOLADJ parameter, setting to T by default.\n")
  my.args$WOLADJ <- T
}

setDefault("POOLED",F)
# If POOLED is "0" or "1" map these to "FALSE" and "TRUE"
if(is.character(my.args$POOLED)) {
  if(my.args$POOLED == "0"){my.args$POOLED <- F}
  if(my.args$POOLED == "1"){my.args$POOLED <- T}
}
# Convert any other type to logical
if(!is.logical(my.args$POOLED)){my.args$POOLED <- as.logical(my.args$POOLED)}
# If conversion produced NA, just set to T
if(is.na(my.args$POOLED)) {
  cat("WARNING: Did not recognize value of POOLED parameter, setting to T by default.\n")
  my.args$POOLED <- T
}

setDefault("TRMT","")

setDefault("VRFILE","")

setDefault("PERM",as.integer(NA))
if(!is.na(my.args$PERM)) {
  # If PERM is character, convert to integer
  if(is.character(my.args$PERM)){my.args$PERM <- as.integer(my.args$PERM)}
  # Round numeric values
  if(is.numeric(my.args$PERM)){
    cat("WARNING: Rounding PERM=", my.args$PERM, " to nearest integer.\n", sep="")
    my.args$PERM <- as.integer(round(my.args$PERM))
  }
  # Don't allow negative values here!
  if(my.args$PERM < 0){
    stop("Cannot have negative PERM index.")
  }
  
  # If PERM specified, seed the RNG
  # Use multiples of large numbers, then take modulo to get within integer range
  my.seed <- as.integer(round((8210678 * (my.args$PERM ** 2.5)) %% (10**9)))
  cat("Seeding RNG with", my.seed, "\n")
  set.seed(my.seed)
}

setDefault("CORES",as.integer(NA))
if(is.na(my.args$CORES)) {
  # If NA, see if environment variable SLURM_JOB_CPUS_PER_NODE exists
  slurm.cpus <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE")
  if(slurm.cpus != "") {
    slurm.cpus <- as.integer(slurm.cpus)
    if(slurm.cpus > 0) {
      my.args$CORES <- slurm.cpus
      cat("Inferring CORES =", slurm.cpus, "based on SLURM_JOB_CPUS_PER_NODE.\n")
    }
  }
}
# If STILL NA, use 1 by default
if(is.na(my.args$CORES)) {
  my.args$CORES <- 1
  cat("Setting CORES = 1 by default.\n")
}
# If CORES = ALL, use detectCores()-1
if(is.character(my.args$CORE) & (my.args$CORE == "ALL")) {
  my.args$CORES <- as.integer(detectCores()-1)
  cat("Setting CORES =", my.args$CORES, "\n")
}
# Otherwise, attempt to convert to integer
if(!is.integer(my.args$CORE)) {
  my.args$CORE <- as.integer(my.args$CORE)
  cat("Setting CORES =", my.args$CORE, "\n")
}

my.args$PSEUDO <- pseudo.count
my.args$LINEINFO <- line.file
my.args$FDR <- fdr.cutoff

# Extract path to EXPR input file, use it as the output dir
my.args$OUTDIR <- paste0(dirname(my.args$EXPR), "/")
# If PERM is specified, append "Perm" and make sure this subdir exists
if(!is.na(my.args$PERM)) {
  my.args$OUTDIR <- paste0(my.args$OUTDIR, "Perm/")
  if(!file.exists(my.args$OUTDIR)) {
    dir.create(path=my.args$OUTDIR, recursive=T, showWarnings = F)
  }
}

# Base name of COUNTS file comes from shaving off subdirectory part of path, and trailing ".txt"
my.args$BASENAME <- sub("[.]txt$", "", basename(my.args$EXPR))

# OUTSTUB will be OUTDIR/BASENAME(_NQ)(_Wol)
# This is used to determine output file paths
my.args$OUTSTUB <- paste0(my.args$OUTDIR, my.args$BASENAME)
if(my.args$NQ) {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_NQ")
}
if(my.args$VRFILE != "") {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_VR")
}
# Wol is not longer flag in outstub, it's only used in the line means file 
# if(my.args$WOLADJ) {
#  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_Wol")
# } else {
#  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_noWol")
# }
if(!is.na(my.args$PERM)) {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_Perm", my.args$PERM)
}

# H2 output table will be OUTSTUB_H2.txt
# Just overwrite if they exist - these output files are tracked by git so it's fine
my.args$H2FILE <- my.args$OUTSTUB
if(my.args$POOLED) {
  my.args$H2FILE <- paste0(my.args$H2FILE, "_Pooled")
}
my.args$H2FILE <- paste0(my.args$H2FILE, "_model_results.txt")

# Suffix of BASENAME is used to infer type of data, but is left in the basename so that this script can be tested on several different normalization strategies
if(grepl("[_]counts$", my.args$BASENAME)) {
  my.args$TYPE <- "counts"
  # my.args$BASENAME <- sub("[_]counts$", "", my.args$BASENAME)
  cat("WARNING: Input file appears to be raw count data, will transform to log2 scale, but this script is intended for normalized expression (rpm, fpkm, or rle)!\n")
} else if(grepl("[_]rpm$", my.args$BASENAME)) {
  my.args$TYPE <- "rpm"
  # my.args$BASENAME <- sub("[_]rpm$", "", my.args$BASENAME)
  cat("Treating input expression file as reads per million (rpm), will transform to log2 scale before fitting lme.\n")
} else if(grepl("[_]fpkm$", my.args$BASENAME)) {
  my.args$TYPE <- "fpkm"
  # my.args$BASENAME <- sub("[_]fpkm$", "", my.args$BASENAME)
  cat("Treating input expression file as features per kb per million (fpkm), will transform to log2 scale before fitting lme.\n")
} else if(grepl("[_]rle$", my.args$BASENAME)) {
  my.args$TYPE <- "rle"
  # my.args$BASENAME <- sub("[_]rle$", "", my.args$BASENAME)
  cat("Treating input expression file as relative log expression (rle), will NOT apply filtering or log2 transformation.\n")
} else {
  my.args$TYPE <- "unknown"
  cat("WARNING: Could not infer expression file type based on filename suffix, will transform to log2 scale before fitting lme.\n")
}
cat("\n")

# Set up parallel backend with desired number of CPUs
cat("Using", my.args$CORES, "CPUs for parallel tasks.\n")
registerDoMC(my.args$CORES)
cat("\n")


# --- LOAD INPUT --- #

# Load the count table
# TO DO: Should encapsulate this in a standard function that loads expr data (counts or normalized) with the optional flag file
expr.table <- read.table(my.args$EXPR, header=T, sep="\t", row.names=1)
# If first column == FLAG, use that to filter the table, then drop the column
dropped.features <- 0
if(colnames(expr.table)[1] == "FLAG") {
  dropped.features <- sum(!(expr.table$FLAG %in% c("OK","LOW","RARE")))
  expr.table <- expr.table[expr.table$FLAG %in% c("OK","LOW","RARE"),]
  expr.table <- expr.table[,2:ncol(expr.table)]
}
# Remove 'X' that R automatically puts in front of column names
colnames(expr.table) <- sub("^X", "", colnames(expr.table))
features <- row.names(expr.table)
cat("Loaded counts for",length(features),"features across",ncol(expr.table),"samples from",my.args$EXPR,"\n")
if(dropped.features > 0) {
  cat("(Dropped",dropped.features,"features based on FLAG column)\n")
}

# Load the master sample table
sample.table <- read.table(my.args$SAMPLES, header=T, sep="\t", row.names=1)
cat("Loaded info on",nrow(sample.table),"samples from",my.args$SAMPLES,"\n")

# How many samples are in common between the two?
all.samples <- intersect(row.names(sample.table),colnames(expr.table))
cat(length(all.samples),"samples matched up between sample info and count table.\n")
if(length(all.samples) < nrow(sample.table)) {
  cat("WARNING: Dropping",length(setdiff(row.names(sample.table),all.samples)),"rows from sample info table.\n")
}
if(length(all.samples) < ncol(expr.table)) {
  cat("WARNING: Dropping",length(setdiff(colnames(expr.table),all.samples)),"columns from count table.\n")
}
sample.table <- sample.table[all.samples,]
expr.table <- expr.table[,all.samples]

# Load the LINE data
line.table <- read.table(my.args$LINEINFO, header=T, sep="\t", row.names=1)
cat("Loaded info on",nrow(line.table),"lines from",my.args$LINEINFO,"\n")
# Remove DGRP_ from row names
row.names(line.table) <- sub("^DGRP[_]","",row.names(line.table))
# Make sure all lines in sample.table are in line.table
stopifnot(all(as.character(sample.table$LINE) %in% row.names(line.table)))
cat(sum(line.table$Wolbachia),"lines have Wolbachia.\n")
sample.table$WOLBACHIA <- line.table[as.character(sample.table$LINE),"Wolbachia"]
cat(round(mean(sample.table$WOLBACHIA)*100),"% of samples have Wolbachia.\n")
cat("\n")

# ADDED 9/20/16 - optional Treatment term
# Determine whether there is Treatment column
use.trmt <- F
if(my.args$TRMT == "") {
  # Add a TRMT column with 1 in all values
  sample.table[,"TRMT"] <- 1
} else {
  # Confirm that the TRMT column is present in sample.table
  if(my.args$TRMT %in% colnames(sample.table)) {
    use.trmt <- T
    cat("Using", my.args$TRMT, "column for modeling treatment effects and interactions.\n")
    # Treatment column should have exactly two groups:
    trmt.groups <- unique(sample.table[,my.args$TRMT])
    if(length(trmt.groups) != 2) {
      stop(my.args$TRMT, " contains ", length(trmt.groups), " groups, this script currently only handles two treatment groups.\n")
    }
    # If values are not 0/1 integer, then re-map here
    if(is.integer(sample.table[,my.args$TRMT])) {
      cat("0 = Control\t1 = Treated\n")
      # Just copy the column to TRMT as is
      sample.table[,"TRMT"] <- sample.table[,my.args$TRMT]
    } else if(is.logical(sample.table[,my.args$TRMT])) {
      cat("FALSE = Control\tTRUE = Treated\n")
      sample.table[,"TRMT"] <- as.integer(sample.table[,my.args$TRMT])
    } else {
      # Any other value type, map the first value to 0, second value to 1
      cat(trmt.groups[1], "= Control\t", trmt.groups[2], "= Treated\n")
      sample.table[,"TRMT"] <- as.integer(factor(sample.table[,my.args$TRMT], levels=trmt.groups))-1
    }
    cat("\n")
  } else {
    stop("Specified Treatment column ", my.args$TRMT, " was not found in ", my.args$SAMPLES, "\n")
  }
}

# Filter rarely expressed features, convert to log2 if not already rle data...
if(my.args$TYPE != "rle") {
  min.rep.expr.perc <- apply(expr.table, 1, function(x){mean(x > 0)})
  
  non.rep.features <- names(min.rep.expr.perc)[min.rep.expr.perc < 0.05]
  
  cat("Removing", length(non.rep.features), "features with expression in <",my.args$FDR*100,"% of all samples.\n")
  keep.features <- setdiff(row.names(expr.table), non.rep.features)
  expr.table <- expr.table[keep.features,]
  
  cat("Converting to log2 scale, with pseudo-count =", my.args$PSEUDO, "\n")
  expr.table <- log2(expr.table + my.args$PSEUDO)
} else {
  cat("No filtering or log2 transformation applied, data is assumed to be log2 scale already (rle type).\n")
}

# If specified, load the variant rate file
vr.table <- NULL
if(my.args$VRFILE != "") {
  stopifnot(file.exists(my.args$VRFILE))
  cat("Loading variant rate data from:", my.args$VRFILE, "\n")
  vr.table <- read.table(my.args$VRFILE, header=T, sep="\t", row.names=1)
  cat("Loaded variant rates for", nrow(vr.table), "genes in", length(setdiff(colnames(vr.table), "TOTAL")), "lines.\n")
  keep.genes <- row.names(expr.table) %in% row.names(vr.table)
  if(!all(keep.genes)) {
    cat("WARNING: Missing variant rates for", sum(!keep.genes), "genes.\n")
    cat("These genes will be dropped and H2 model will only run for", sum(keep.genes), "genes.\n")
    expr.table <- expr.table[keep.genes,]
  }
  # Drop variant rates for genes missing from expr table (don't need to warn on this one)
  stopifnot(all(row.names(expr.table) %in% row.names(vr.table)))
  vr.table <- vr.table[row.names(expr.table),]
}

# Check if this is combined known/novel gene table:
split.known.novel <- F
if(any(grepl("^FBgn",row.names(expr.table))) & any(grepl("^XLOC",row.names(expr.table)))) {
  split.known.novel <- T
  known.rows <- grep("^FBgn",row.names(expr.table),value=T)
  novel.rows <- grep("^XLOC",row.names(expr.table),value=T)
}

# TEMP DEBUGGING
# expr.table <- expr.table[c(known.rows[1:500], novel.rows[1:500]),]


# NOTE: I'm assuming there will be huge sex effects in the microbiome and transposon results, but that may not be true
# assess_heritability_pooled.R is used to test that assumption and compute heritability and line effects across sexes

# If SEX column exists in sample.table, split the expr table by sex
# Also define the groups that will be used for line means
sample.subsets <- list()
if(my.args$POOLED) {
  if("SEX" %in% colnames(sample.table)) {
    all.lines <- unique(sample.table$LINE)
    female.line.reps <- unlist(lapply(all.lines, function(x){sum((sample.table$LINE==x) & (sample.table$SEX=="F"))}))
    names(female.line.reps) <- all.lines
    male.line.reps <- unlist(lapply(all.lines, function(x){sum((sample.table$LINE==x) & (sample.table$SEX=="M"))}))
    names(male.line.reps) <- all.lines
    complete.lines <- all.lines[(female.line.reps >= 2) & (male.line.reps >= 2)]
    drop.lines <- setdiff(all.lines, complete.lines)
    cat("Dropping lines with incomplete 2F+2M pairing:", drop.lines, "\n")
    keep.samples <- row.names(sample.table)[sample.table$LINE %in% complete.lines]
    cat("Fitting pooled sex by line model using", length(keep.samples), "covering", length(complete.lines), "lines.\n")
    sample.subsets[["Pooled"]] <- keep.samples
    if(use.trmt) {
      sample.table$GROUP <- paste(sample.table$LINE, sample.table$SEX, sample.table[,my.args$TRMT], sep="_")
    } else {
      sample.table$GROUP <- paste(sample.table$LINE, sample.table$SEX, sep="_")
    }
  } else {
    stop("No SEX column in", my.args$SAMPLES, "- cannot run pooled sex by line models.\n")
  }
} else {
  if("SEX" %in% colnames(sample.table)) {
    cat("Splitting expression samples by sex.\n")
    # TO DO: SHOULD DROP LINES WITH SINGLE REP WITHIN SEX HERE
    # SHOULD BE ABLE TO DO BY IDENTIFYING GROUPS WITH SINGLE REP WITHIN EACH SEX?
    for(sex in sort(unique(sample.table$SEX))) {
      sample.subsets[[sex]] <- row.names(sample.table)[sample.table$SEX==sex]
    }
    if(use.trmt) {
      sample.table$GROUP <- paste(sample.table$LINE, sample.table$SEX, sample.table[,my.args$TRMT], sep="_")
    } else {
      sample.table$GROUP <- paste(sample.table$LINE, sample.table$SEX, sep="_")
    }
  } else {
    cat("No SEX column in", my.args$SAMPLES,"- will compute single sex linear models across all samples.\n")
    sample.subsets[["All"]] <- row.names(sample.table)
    # Add a SEX column with constant value (U or Unknown)
    sample.table[,"SEX"] <- "U"
    if(use.trmt) {
      sample.table$GROUP <- as.character(sample.table$LINE, sample.table[,my.args$TRMT], sep="_")
    } else {
      sample.table$GROUP <- as.character(sample.table$LINE)
    }
  }
}


# --- Core Sub-Routines --- #

# NOTE: Splitting out the individual model-fitting steps into sub-routines seemed to 
# cause problems because passing a complete model result back from a function appears to lose some critical data
# This is a flaw in the lmer module design as best I can tell

# This function combines one expression vector with all model data
# i = the gene index to extract
# my.expr should be a matrix with columns = genes, samples = rows (should also be scaled and normalized already)
# my.samples should be a data frame with samples = rows, and columns for SEX, WOLBACHIA, LINE, TRMT, and GROUP
# When there is not TRMT or SEX term, these should be present but set to a fixed value
# GROUP is not used in models, but is used to label line means
getModelInput <- function(i, my.expr, my.samples, var.rate=NULL) {
  stopifnot(nrow(my.expr)==nrow(my.samples))
  stopifnot(row.names(my.expr)==row.names(my.samples))
  my.model.input <- cbind(EXPR=my.expr[,i], my.samples[,c("SEX","WOLBACHIA","LINE","TRMT","GROUP")])
  if(!is.null(var.rate)) {
    # Extract variant rate data if provided
    stopifnot(all(dim(my.expr)==dim(var.rate)))
    stopifnot(all(row.names(my.expr)==row.names(var.rate)))
    stopifnot(all(colnames(my.expr)==colnames(var.rate)))
    my.model.input[,"VR"] <- var.rate[,i]
  }
  return(my.model.input)
}

# This function takes a single model input data frame and performs ALL modeling steps
# model.input should be a data frame as returned by getModelInput above
# Set verbose=T to see progress reporting - DO NOT DO FOR FULL/PARALLEL PROCESSING!
# All other params currently pulled from global my.args
# Return type is a list with two members:
# $model.results = vector of all parameters for model result table (except adjusted p-values)
# $line.means = vector of line means
# NOTE: This function does not handle SCALING any model params or line means back to the original distributions
#       The model.input is taken at face value, without any scaling data
#       This function also does not handle p-value adjustments, since it only sees one gene expression vector at a time
# TO DO: Could restructure this to fit all models at the beginning (there's probably some redundancy in current code)
#     Then, depending on which params are set, can decide which tests anova and rand tests to do?
# TO DO: Could set all parameters explicitly, rather than relying on globals
#     (only a few are actually used here: my.args$POOLED, my.args$WOLADJ, use.trmt)
modelGeneticVar <- function(model.input, verbose=F) {
  # This vector will store named model result values:
  my.results <- c()
  
  # If there's a VR column present, regress out these effects first
  if("VR" %in% colnames(model.input)) {
    # Regress out by fitting a standard LM first
    # Residuals will be the new expr for main model
    # NOTE: We *COULD* try incorporating these directly in the main model, but I think this is the better way to do things
    # This pulls out the maximum amount of variance attributable to alignment bias BEFORE fitting the main model of interest
    if(verbose) cat("Correcting for gene level variation rate.")
    # TO DO: For POOLED or TRMT model, could include ixn terms here?
    #   This would be for genes that are only expressed in one sex or one treatment group
    #   But I don't think this is necessary - we don't really care about this correction in either of those contexts yet
    #   Would also need to add back the primary SEX and/or TRMT effects here
    vr.lm <- lm(EXPR ~ VR, data=model.input)
    
    # Compute Coefficient and Avg Effect Size (This can be scaled back later)
    my.results["VR.Coeff"] <- vr.lm$coefficients["VR"]
    my.results["VR.Log2FC"] <- mean(vr.lm$coefficients["VR"] * model.input$VR)
    # These can be NA if this gene has no recorded variants, so VR is 0 for all
    if(is.na(my.results["VR.Coeff"])) my.results["VR.Coeff"] <- 0
    if(is.na(my.results["VR.Log2FC"])) my.results["VR.Log2FC"] <- 0
    
    # Compute the corrected expression values
    # Add intercept back into each column of resid to keep similar scale
    adj.expr <- resid(vr.lm) + vr.lm$coefficients["(Intercept)"]
    
    orig.var <- var(model.input$EXPR)
    adj.var <- var(adj.expr)
    my.results["VR.Perc.Var"] <- (orig.var - adj.var) / orig.var
    
    # Use ANOVA to get P-value
    my.results["VR.PVal"] <- anova(vr.lm)$Pr[1]
    # If NA, set to P-value=1
    if(is.na(my.results["VR.PVal"])) my.results["VR.PVal"] <- 1
    
    # Replace the original expression data with the adjusted expression data
    model.input$EXPR <- adj.expr
  }

  # Fit model w/out Wolbachia (within-sex or pooled depending on command-line params)
  # Which model depends on whether we're doing Pooled-Sex and/or using a Treatment term
  if(my.args$POOLED) {
    if(verbose) cat("Fitting Pooled-Sex Mixed Effects Model without Wolbachia adjustment...\n")
    if(use.trmt) {
      my.noWol.lmer <- lmer(EXPR ~ SEX + TRMT + (1|LINE) + (1|SEX:LINE) + (1|TRMT:LINE) + (1|SEX:TRMT:LINE), data=model.input)
      # This might fail - this might be TRMT(second value), or their could be multiple if there 3+ TRMT groups...
      my.results["NoWol.Trmt.Log2FC"] <- fixef(my.noWol.lmer)["TRMT"]
    } else {
      my.noWol.lmer <- lmer(EXPR ~ SEX + (1|LINE) + (1|SEX:LINE), data=model.input)
    }
    my.results["NoWol.Sex.Log2FC"] <- fixef(my.noWol.lmer)["SEXM"]
  } else {  
    if(verbose) cat("Fitting Within-Sex Random Line Effects Model without Wolbachia adjustment...\n")
    if(use.trmt) {
      my.noWol.lmer <- lmer(EXPR ~ TRMT + (1|LINE) + (1|TRMT:LINE), data=model.input)
      my.results["NoWol.Trmt.Log2FC"] <- fixef(my.noWol.lmer)["TRMT"]
    } else {
      my.noWol.lmer <- lmer(EXPR ~ (1|LINE), data=model.input)
    }
  }
  # Extract variance components and compute H2
  # NOTE: No longer scaling back to the original expression variance, but this cancels out in H2 calculation anyway
  my.noWol.varcorr <- as.data.frame(VarCorr(my.noWol.lmer))
  my.noWol.varcorr.names <- my.noWol.varcorr$grp
  my.noWol.varcorr <- my.noWol.varcorr$vcov
  names(my.noWol.varcorr) <- my.noWol.varcorr.names
  # genet.var term will sum up ALL genetic variance components
  my.noWol.genet.var <- my.noWol.varcorr["LINE"]
  for(opt.term in c("SEX:LINE","TRMT:LINE","SEX:TRMT:LINE")) {
    if(opt.term %in% names(my.noWol.varcorr)) {
      my.noWol.genet.var <- my.noWol.genet.var + my.noWol.varcorr[opt.term]
    }
  }
  # env.var is just the residual variance
  my.noWol.env.var <- my.noWol.varcorr["Residual"]
  my.results["NoWol.H2"] <- my.noWol.genet.var / (my.noWol.genet.var + my.noWol.env.var)
  # Test the random effect using lmerTest package and extract p-values to my.results vector
  my.noWol.rand <- rand(my.noWol.lmer)$rand.table
  my.results["NoWol.Line.PVal"] <- my.noWol.rand["LINE","p.value"]
  if("SEX:LINE" %in% row.names(my.noWol.rand)) {
    my.results["NoWol.Sex.Line.PVal"] <- my.noWol.rand["SEX:LINE","p.value"]
  }
  if("TRMT:LINE" %in% row.names(my.noWol.rand)) {
    my.results["NoWol.Trmt.Line.PVal"] <- my.noWol.rand["TRMT:LINE","p.value"]
  }
  if("SEX:TRMT:LINE" %in% row.names(my.noWol.rand)) {
    my.results["NoWol.Sex.Trmt.Line.PVal"] <- my.noWol.rand["SEX:TRMT:LINE","p.value"]
  }
  
  # Now fit model with Wolbachia correction
  if(my.args$POOLED) {
    if(verbose) cat("Fitting Pooled-Sex Mixed Effects Model with Wolbachia corrections...\n")
    if(use.trmt) {
      my.adjWol.lmer <- lmer(EXPR ~ SEX + WOLBACHIA + SEX:WOLBACHIA + TRMT + (1|LINE) + (1|SEX:LINE) + (1|TRMT:LINE) + (1|SEX:TRMT:LINE), data=model.input)
      my.results["WolAdj.Trmt.Log2FC"] <- fixef(my.adjWol.lmer)["TRMT"]
    } else {
      my.adjWol.lmer <- lmer(EXPR ~ SEX + WOLBACHIA + SEX:WOLBACHIA + (1|LINE) + (1|SEX:LINE), data=model.input)
    }
    # OLD METHOD FOR GETTING SEX EFFECT P-VALUE
    # # Also fit model without sex effects or interactions for anova p-values below
    # # This can be fit with ML method from the start, since it is only used for ANOVA
    # if(verbose) cat("Fitting Pooled-Sex Mixed Effects Model with Wolbachia corrections but NO sex effects...\n")
    # if(use.trmt) {
    #  my.adjWol.noSex.lmer.ml <- lmer(EXPR ~ WOLBACHIA + TRMT + (1|LINE) + (1|TRMT:LINE), data=model.input, REML=F)
    # } else {
    #  my.adjWol.noSex.lmer.ml <- lmer(EXPR ~ WOLBACHIA + (1|LINE), data=model.input, REML=F)
    # }
    # if(verbose) cat("Comparing models with and without Sex effects by ANOVA...\n")
    # # Manually refitting by ML prevents extensive warning messages
    # # my.adjWol.noSex.lmer.ml <- refitML(my.adjWol.noSex.lmer)
    # my.adjWol.lmer.ml <- refitML(my.adjWol.lmer)
    # my.results["WolAdj.Sex.PVal"] <- anova(my.adjWol.noSex.lmer.ml, my.adjWol.lmer.ml)$Pr[2]
    # NEW METHOD - JUST USE ANOVA FUNCTION
    my.adjWol.anova <- anova(my.adjWol.lmer)
    if("Pr(>F)" %in% colnames(my.adjWol.anova)) {
      my.results["WolAdj.Sex.PVal"] <- my.adjWol.anova["SEX","Pr(>F)"]
    } else {
      my.results["WolAdj.Sex.PVal"] <- NA
    }
    
    # When treatment term included, try fitting without that as well for ANOVA test
    if(use.trmt) {
      # OLD METHOD - SKIP
      # if(verbose) cat("Fitting Pooled-Sex Mixed Effects Model with Wolbachia corrections but NO treatment effects...\n")
      # my.adjWol.noTrmt.lmer.ml <- lmer(EXPR ~ SEX + WOLBACHIA + SEX:WOLBACHIA + (1|LINE) + (1|SEX:LINE), data=model.input, REML=F)
      # # TO DO: Not sure if there are other values in this probability vector that should be used instead?
      # my.results["WolAdj.Trmt.PVal"] <- anova(my.adjWol.noTrmt.lmer.ml, my.adjWol.lmer.ml)$Pr[2]
      if("Pr(>F)" %in% colnames(my.adjWol.anova)) {
        my.results["WolAdj.Trmt.PVal"] <- my.adjWol.anova["TRMT","Pr(>F)"]
      } else {
        my.results["WolAdj.Trmt.PVal"] <- NA
      }
    }
  } else {
    if(verbose) cat("Fitting Within-Sex Mixed Effects Model with Fixed Wolbachia Effects and Random Line Effects...\n")
    if(use.trmt) {
      my.adjWol.lmer <- lmer(EXPR ~ WOLBACHIA + TRMT + (1|LINE) + (1|TRMT:LINE), data=model.input)
      my.results["WolAdj.Trmt.Log2FC"] <- fixef(my.adjWol.lmer)["TRMT"]
      
      # OLD APPROACH - SKIP
      # NOTE: THE CORRECT WAY TO DO THIS WOULD HAVE BEEN TO FIT MODEL THAT IS ONLY MISSING TRMT (BUT STILL HAS TRMT:LINE)
      # # When treatment term included, try fitting without that as well for ANOVA test
      # if(verbose) cat("Fitting Within-Sex Mixed Effects Model with Fixed Wolbachia corrections but NO treatment effects...\n")
      # my.adjWol.noTrmt.lmer.ml <- lmer(EXPR ~ WOLBACHIA + (1|LINE), data=model.input, REML=F)
      # # Refit the original treatment model with ML method, which is better for ANOVA test
      # my.adjWol.lmer.ml <- refitML(my.adjWol.lmer)
      # # TO DO: Not sure if there are other values in this probability vector that should be used instead?
      # my.results["WolAdj.Trmt.PVal"] <- anova(my.adjWol.noTrmt.lmer.ml, my.adjWol.lmer.ml)$Pr[2]
      
      # NEW APPROACH
      # This is a true F-Test although the results are not much different?
      # NOTE: The first entry in this vector is the Wolbachia fixed effect p-value
      # TO DO: Drop all the fitting of Wolbachia+/- model, just use this test
      my.results["WolAdj.Trmt.PVal"] <- anova(my.adjWol.lmer)["TRMT","Pr(>F)"]
    } else {
      my.adjWol.lmer <- lmer(EXPR ~ WOLBACHIA + (1|LINE), data=model.input)
    }
  }
  # Extract the fixed effects
  my.adjWol.fixef <- fixef(my.adjWol.lmer)
  my.adjWol.wolEff <- my.adjWol.fixef["WOLBACHIA"]
  if("SEXM" %in% names(my.adjWol.fixef)) {
    my.results["WolAdj.Sex.Log2FC"] <- my.adjWol.fixef["SEXM"]
  }
  if("TRMT" %in% names(my.adjWol.fixef)) {
    my.results["WolAdj.Trmt.Log2FC"] <- my.adjWol.fixef["TRMT"]
  }
  
  # Extract the variance components from the full model:
  my.adjWol.varcorr <- as.data.frame(VarCorr(my.adjWol.lmer))
  my.adjWol.varcorr.names <- my.adjWol.varcorr$grp
  my.adjWol.varcorr <- my.adjWol.varcorr$vcov
  names(my.adjWol.varcorr) <- my.adjWol.varcorr.names
  # genet.var term will sum up ALL genetic variance components
  my.adjWol.genet.var <- my.adjWol.varcorr["LINE"]
  for(opt.term in c("SEX:LINE","TRMT:LINE","SEX:TRMT:LINE")) {
    if(opt.term %in% names(my.adjWol.varcorr)) {
      my.adjWol.genet.var <- my.adjWol.genet.var + my.adjWol.varcorr[opt.term]
    }
  }
  # env.var is just the residual variance
  my.adjWol.env.var <- my.adjWol.varcorr["Residual"]
  my.results["WolAdj.H2"] <- my.adjWol.genet.var / (my.adjWol.genet.var + my.adjWol.env.var)
  
  # Extract the random effect p-values using lmerTest package
  my.adjWol.rand <- rand(my.adjWol.lmer)$rand.table
  my.results["WolAdj.Line.PVal"] <- my.adjWol.rand["LINE","p.value"]
  if("SEX:LINE" %in% row.names(my.adjWol.rand)) {
    my.results["WolAdj.Sex.Line.PVal"] <- my.adjWol.rand["SEX:LINE","p.value"]
  }
  if("TRMT:LINE" %in% row.names(my.adjWol.rand)) {
    my.results["WolAdj.Trmt.Line.PVal"] <- my.adjWol.rand["TRMT:LINE","p.value"]
  }
  if("SEX:TRMT:LINE" %in% row.names(my.adjWol.rand)) {
    my.results["WolAdj.Sex.Trmt.Line.PVal"] <- my.adjWol.rand["SEX:TRMT:LINE","p.value"]
  }
  
  # Compare the two models for overall Wolbachia effect:
  if(verbose) cat("Comparing Wolbachia-free and Wolbachia-adjusted models by ANOVA...\n")
  # OLD METHOD - SKIP
  # # Manually refitting by ML to avoid extensive warning messages
  # my.noWol.lmer.ml <- refitML(my.noWol.lmer)
  # my.adjWol.lmer.ml <- refitML(my.adjWol.lmer)
  # my.wol.anova.pval <- anova(my.noWol.lmer.ml, my.adjWol.lmer.ml)$Pr[2]
  # NEW METHOD - Just use standard Type III Anova w/ Satterthwaite approximation of d.f.
  my.wol.anova.pval <- anova(my.adjWol.lmer)["WOLBACHIA","Pr(>F)"]
  if(is.null(my.wol.anova.pval)) {
    my.wol.anova.pval <- NA
  }
  
  # To compute Perc.Var, need to compute the residual values from each model and look at the difference in variance
  # TO DO: FIX THIS: I don't think this is correct - could argue that should just take the original expr variance compared to subtracting the Wolbachia fixed effects
  # DO NOT need to add back the model intercept or completely unscale the values
  # TO DO: MIGHT BE USEFUL TO MAKE THIS A SUB-ROUTINE
  my.noWol.resid.var <- var(resid(my.noWol.lmer))
  my.adjWol.resid.var <- var(resid(my.adjWol.lmer))
  my.wol.perc.var <- (my.noWol.resid.var - my.adjWol.resid.var) / my.noWol.resid.var
  
  # Put these in result table
  my.results["Wol.Perc.Var"] <- my.wol.perc.var
  my.results["Wol.Log2FC"] <- my.adjWol.wolEff
  my.results["Wol.PVal"] <- my.wol.anova.pval
  
  # Compute line means
  if(verbose) cat("Extracting line means using ")
  # Which model depends on on my.args$WOLADJ, my.args$POOLED, and use.trmt
  # But in all cases, EXCLUDE the Wolbachia effect
  if(my.args$POOLED) {
    if(my.args$WOLADJ) {
      if(verbose) cat("Pooled-Sex Wolbachia-adjusted mixed effects model.\n")
      if(use.trmt) {
        my.fitted <- predict(my.adjWol.lmer, re.form=~SEX+TRMT+(1|LINE)+(1|SEX:LINE)+(1|TRMT:LINE)+(1|SEX:TRMT:LINE))
      } else {
        my.fitted <- predict(my.adjWol.lmer, re.form=~SEX+(1|LINE)+(1|SEX:LINE))
      }
    } else {
      if(verbose) cat("Pooled-Sex mixed effects model with NO Wolbachia correction.\n")
      if(use.trmt) {
        my.fitted <- predict(my.noWol.lmer, re.form=~SEX+TRMT+(1|LINE)+(1|SEX:LINE)+(1|TRMT:LINE)+(1|SEX:TRMT:LINE))
      } else {
        my.fitted <- predict(my.noWol.lmer, re.form=~SEX+(1|LINE)+(1|SEX:LINE))
      }
    }
  } else {
    if(my.args$WOLADJ) {
      if(verbose) cat("Within-Sex Wolbachia-adjusted mixed effects model.\n")
      if(use.trmt) {
        my.fitted <- predict(my.adjWol.lmer, re.form=~TRMT+(1|LINE)+(1|TRMT:LINE))
      } else {
        my.fitted <- predict(my.adjWol.lmer, re.form=~(1|LINE))
      }
    } else {
      if(verbose) cat("Within-Sex random effects model with NO Wolbachia correction.\n")
      if(use.trmt) {
        my.fitted <- predict(my.noWol.lmer, re.form=~TRMT+(1|LINE)+(1|TRMT:LINE))
      } else {
        my.fitted <- predict(my.noWol.lmer, re.form=~(1|LINE))
      }
    }
  }
  
  # Create a list structure that maps Groups to Samples, with groups ordered by line ID
  my.groups <- unique(model.input$GROUP[order(model.input$TRMT, model.input$SEX, model.input$LINE, decreasing=F)])
  names(my.groups) <- my.groups
  my.groups <- lapply(my.groups, function(x){row.names(model.input)[model.input$GROUP==x]})
  
  # Loop over groups, make sure all fitted values match, then collapse into table
  my.line.means <- unlist(lapply(my.groups, function(x){
    if(length(x) == 1) {
      return(my.fitted[x])
    } else {
      stopifnot(all(my.fitted[x] == my.fitted[x[1]]))
      return(mean(my.fitted[x]))
    }
  }))
  names(my.line.means) <- names(my.groups)
  
  # RETURN LIST WITH MEMBERS model.results, line.means
  return(list(model.results=my.results, line.means=my.line.means))
}

# Function that takes a list of results from modelGeneticVar and compiles the model result or line mean table
# List member names should be the Gene/Feature IDs
# Set extract=line.means to get line means table, default is model.results table
extractModelTable <- function(model.outputs, extract="model.results") {
  my.table <- foreach(i=1:length(model.outputs), .combine='rbind') %do% {model.outputs[[i]][[extract]]}
  if(!is.null(names(model.outputs))) {
    row.names(my.table) <- names(model.outputs)
  }
  return(as.data.frame(my.table))
}

# Function that inserts a multiple-testing corrected column immediately after P-value column
correctPvalCols <- function(
  my.table,       # data.frame with p-values in pval.cols
  pval.cols=c(),  # must be by name, if empty default to all columns ending in [.]PVal 
  method="BH",    # Method for correcting the p-values
  suffix=""       # Suffix to add to adjusted p-value column names (leave blank for auto rules)
) {
  stopifnot(!is.null(colnames(my.table)))
  if(length(pval.cols)==0) {
    pval.cols <- grep("[.]PVal", colnames(my.table), value=T)
  }
  stopifnot(length(pval.cols) > 0)
  stopifnot(all(pval.cols %in% colnames(my.table)))
  if(suffix=="") {
    if(method=="bonferroni") {
      suffix <- ".BF"
    } else if(method %in% c("fdr","BH")) {
      suffix <- ".FDR"
    } else {
      suffix <- paste0(".", method)
    }
  }
  
  # Loop over each pval.col
  for(j in pval.cols) {
    # Map name to col number
    j <- which(colnames(my.table)==j)[1]
    # Do multiple testing correction
    j.adj <- p.adjust(my.table[,j], method=method)
    # Insert the corrected column immediately after j
    if(j < ncol(my.table)) {
      my.table <- cbind(my.table[,1:j], j.adj, my.table[,(j+1):ncol(my.table)])
    } else {
      my.table <- cbind(my.table, j.adj)
    }
    colnames(my.table)[j+1] <- paste0(colnames(my.table)[j], suffix)
  }
  
  # Return the expanded table
  return(my.table)
}



# --- Main Loop --- #

# Loop over sample subsets
cat("Running lme4 model-fitting on", length(sample.subsets), "sample subsets.\n\n")
ss.results <- list()
for(ss in names(sample.subsets)) {
  
  cat("- Running lme4 model-fitting on",ss,"samples ")
  if(my.args$NQ) {
    cat("w/ Normal Quantile Transform ")
  }
  cat("-\n\n")
  
  # Extract the subset-specific expression
  ss.expr <- expr.table[,sample.subsets[[ss]]]
  ss.samples <- sample.table[sample.subsets[[ss]],]
  
  # Drop features with 0 variance
  drop.features <- (apply(ss.expr, 1, var) < (10**-8))
  cat("Dropping", sum(drop.features), "features with no variance in this sample subset.\n")
  ss.expr <- ss.expr[!drop.features,]
  
  # NON-RLE ONLY
  if(my.args$TYPE != "rle") {
    min.rep.expr.perc <- apply(ss.expr, 1, function(x){mean(x > log2(my.args$PSEUDO))})
    non.rep.features <- names(min.rep.expr.perc)[min.rep.expr.perc < 0.05]
    
    cat("Removing", length(non.rep.features), "features with expression in <",my.args$FDR*100,"% of this subset of samples.\n")
    keep.features <- setdiff(row.names(ss.expr), non.rep.features)
    ss.expr <- ss.expr[keep.features,]
  }
  cat("\n")
  
  # NEW - 7/8/16
  # When PERM optin specified, 
  # Permute the sample IDs to get bg distribution of H2 values and p-values
  # TO DO: THIS IS NOT VALID FOR POOLED MODEL!!!!
  #   A better way to do this would be to define a bunch of possible permutations in a separate script and output to a table
  #   Then load that table here and apply a single permutation 
  if(!is.na(my.args$PERM)) {
    if(my.args$POOLED) {
      stop("PERM method is not yet set up for POOLED model!")
    }
    # Create a random ordering
    cat("Shuffling sample IDs...\n")
    samp.shuf <- sample.int(n=ncol(ss.expr), replace=F)
    real.groups <- lapply(unique(ss.samples$LINE), function(x){which(ss.samples$LINE==x)})
    colnames(ss.expr) <- colnames(ss.expr)[samp.shuf]
    ss.samples <- ss.samples[colnames(ss.expr),]
    perm.groups <- lapply(unique(ss.samples$LINE), function(x){which(ss.samples$LINE==x)})
    # Compute how many biological replicates are still paired together after shuffling
    # Should be 0...
    unperm.groups <- sum(unlist(lapply(perm.groups, function(x){
      max(unlist(lapply(real.groups, function(y){length(intersect(x,y))}))) - 1
    })))
    cat(unperm.groups, "true replicate pairs are still paired after shuffling.\n\n")
  }
  
  # Now transpose the ss.expr data so that features are columns (makes it easier to combine with LINE column for lmer model)
  ss.expr <- t(ss.expr)
  
  if(my.args$NQ) {
    # Normal Quantile Transformation
    cat("Applying Normal Quantile Transformation.\n")
    ss.expr <- apply(ss.expr, 2, function(x){
      qnorm(rank(x, ties.method="average")/(length(x)+1))
    })
  }
  
  # Apply scale function to each expression profile here
  ss.expr <- scale(ss.expr)
  ss.center <- attr(ss.expr, "scaled:center")
  ss.sd <- attr(ss.expr, "scaled:scale")
  # When adding model terms, make sure these terms get added back in at relevant places below!
  
  # TO DO: (LOW) Could put this in a separate script that just does the correction and outputs a new expr table?
  # If vr.table is present, do simple linear regression on that first, and take the residuals
  if(is.null(vr.table)) {
    ss.vr <- NULL
  } else {
    # Just transform vr.table to have same dimensions as ss.expr here
    # Actual modeling of effects will be handled in core modeling function

    # Construct a vector that maps each sample to the appropriate column in vr.table
    ss.samp.vr.col <- paste0("line_", ss.samples$LINE)
    names(ss.samp.vr.col) <- row.names(ss.samples)
    if(!all(ss.samp.vr.col %in% colnames(vr.table))) {
      stop("Missing variant rate data for one or more sample lines.")
    }
    stopifnot(all(names(ss.samp.vr.col)==row.names(ss.expr)))
    # Create variant rate table with same dimensions as ss.expr
    #  - just duplicate line values for replicate samples
    ss.vr <- vr.table[colnames(ss.expr),ss.samp.vr.col]
    colnames(ss.vr) <- names(ss.samp.vr.col)
    ss.vr <- t(ss.vr)
    stopifnot(all(dim(ss.expr)==dim(ss.vr)))
    stopifnot(all(row.names(ss.expr)==row.names(ss.vr)))
    stopifnot(all(colnames(ss.expr)==colnames(ss.vr)))
  }
  
  # NOTE: I tried to modularize the modeling step, but when it's called from within a function
  # it causes problems in some of the downstream functions like rand
  
  # Loop over each gene, run all modeling steps, 
  # and return a list structure with model results and line means
  # This loop can be parallelized and hopefully will be more stable for multi-threading than previous attempts...
  cat("Running core model fitting and analysis on", ncol(ss.expr), "features using", my.args$CORES, "CPUs...\n")
  # Time this step, output CPU usage for this part specifically
  start.time <- proc.time()
  ss.model.output <- foreach(i=1:ncol(ss.expr)) %dopar% {
    i.model.input <- getModelInput(i=i, my.expr=ss.expr, my.samples=ss.samples, var.rate=ss.vr)
    modelGeneticVar(i.model.input)
  }
  end.time <- proc.time()
  names(ss.model.output) <- colnames(ss.expr)
  cat("Model fitting and analysis complete, time usage:\n")
  print(end.time-start.time)
  cat("\n")
  
  # Convert these results to tables
  ss.model.results <- extractModelTable(ss.model.output)
  ss.line.means <- extractModelTable(ss.model.output, extract="line.means")
  
  # Scale certain columns back to original distribution variance/sd
  for(fc.col in c("VR.Coeff", "VR.Log2FC", "Wol.Log2FC", "NoWol.Sex.Log2FC", "WolAdj.Sex.Log2FC", "NoWol.Trmt.Log2FC", "WolAdj.Trmt.Log2FC")) {
    if(fc.col %in% colnames(ss.model.results)) {
      ss.model.results[,fc.col] <- ss.model.results[,fc.col] * ss.sd
    }
  }
  
  # Perform multiple testing correction on all P-value columns
  ss.model.results <- correctPvalCols(ss.model.results, method="BH")
  
  # Report VR correction parameters if computed
  if("VR.Coeff" %in% colnames(ss.model.results)) {
    cat(sum(ss.model.results$VR.PVal.FDR <= my.args$FDR),"features have significant Variation Rate effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
    cat(" =", round(mean(ss.model.results$VR.PVal.FDR <= my.args$FDR)*100), "% of all features tested.\n\n")
    cat("Distribution of % Variance explained by Gene-Level Variation Rate:\n")
    cat("Range:", range(ss.model.results[,"VR.Perc.Var"]), "\n")
    cat("Box Stats:", boxplot.stats(ss.model.results[,"VR.Perc.Var"])$stats,"\n")
    cat("Distribution of Variation Rate Effect Sizes (Log2FC):\n")
    cat("Range:", range(ss.model.results[,"VR.Log2FC"]), "\n")
    cat("Box Stats:", boxplot.stats(ss.model.results[,"VR.Log2FC"])$stats,"\n")
    cat("\n")
  }
  
  # U R HERE - 10/4/16
  # WHY ARE NO FDR VALUES NA - DOES IT TURN THOSE TO 1 AUTOMATICALLY?
  
  # Report how many features are significantly heritable in Wolbachia-free model
  noWol.fdr.cols <- c("NoWol.Line.PVal.FDR", "NoWol.Sex.Line.PVal.FDR", "NoWol.Trmt.Line.PVal.FDR", "NoWol.Sex.Trmt.Line.PVal.FDR")
  noWol.fdr.cols <- intersect(noWol.fdr.cols, colnames(ss.model.results))
  stopifnot(length(noWol.fdr.cols) > 0)
  if(length(noWol.fdr.cols) == 1) {
    signif.rows <- ss.model.results[,noWol.fdr.cols] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    signif.rows[failed.rows] <- F 
  } else {
    signif.rows <- apply(ss.model.results[,noWol.fdr.cols], 1, function(x){any(x <= my.args$FDR, na.rm=F)})
    failed.rows <- apply(ss.model.results[,noWol.fdr.cols], 1, function(x){any(is.na(x))})
  }
  if(sum(failed.rows) > 0) {
    cat(sum(failed.rows),"features failed to converge on random effects model WITHOUT Wolbachia for", ss, "samples\n")
    cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
  }
  cat(sum(signif.rows),"features have significant genetic variance in random effects model WITHOUT Wolbachia for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
  cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n")
  cat("Distribution of H^2 values for significantly heritable features:\n")
  cat("Range:", range(ss.model.results[signif.rows,"NoWol.H2"]), "\n")
  cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"NoWol.H2"])$stats,"\n")
  cat("\n")
  
  adjWol.fdr.cols <- c("WolAdj.Line.PVal.FDR", "WolAdj.Sex.Line.PVal.FDR", "WolAdj.Trmt.Line.PVal.FDR", "WolAdj.Sex.Trmt.Line.PVal.FDR")
  adjWol.fdr.cols <- intersect(adjWol.fdr.cols, colnames(ss.model.results))
  stopifnot(length(adjWol.fdr.cols) > 0)
  if(length(adjWol.fdr.cols) == 1) {
    signif.rows <- ss.model.results[,adjWol.fdr.cols] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    signif.rows[failed.rows] <- F
  } else {
    signif.rows <- apply(ss.model.results[,adjWol.fdr.cols], 1, function(x){any(x <= my.args$FDR, na.rm=F)})
    failed.rows <- apply(ss.model.results[,adjWol.fdr.cols], 1, function(x){any(is.na(x))})
  }
  if(sum(failed.rows) > 0) {
    cat(sum(failed.rows),"features failed to converge on random effects model with Wolbachia effects for", ss, "samples\n")
    cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
  }
  cat(sum(signif.rows),"features have significant genetic variance in mixed effects model with Wolbachia effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
  cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n")
  cat("Distribution of H^2 values for significantly heritable features:\n")
  cat("Range:", range(ss.model.results[signif.rows,"WolAdj.H2"]), "\n")
  cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"WolAdj.H2"])$stats,"\n")
  cat("\n")
  
  signif.rows <- ss.model.results$Wol.PVal.FDR <= my.args$FDR
  failed.rows <- is.na(signif.rows)
  signif.rows[failed.rows] <- F
  if(sum(failed.rows) > 0) {
    cat(sum(failed.rows),"features failed to converge when testing Wolbachia fixed effect ANOVA\n")
    cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
  }
  cat(sum(signif.rows),"features have significant Wolbachia effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
  cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n\n")
  cat("Distribution of Wolbachia % Variance Explained:\n")
  cat("Range:", range(ss.model.results[,"Wol.Perc.Var"]), "\n")
  cat("Box Stats:", boxplot.stats(ss.model.results[,"Wol.Perc.Var"])$stats,"\n")
  cat("Distribution of Wolbachia Effects:\n")
  cat("Range:", range(ss.model.results[,"Wol.Log2FC"]), "\n")
  cat("Box Stats:", boxplot.stats(ss.model.results[,"Wol.Log2FC"])$stats,"\n")
  cat("\n")
  
  if(use.trmt) {
    fdr.col <- "Trmt.PVal.FDR"
    fc.col <- "Trmt.Log2FC"
    if(my.args$WOLADJ) {
      fdr.col <- paste0("WolAdj.", fdr.col)
      fc.col <- paste0("WolAdj.", fc.col)
    } else {
      fdr.col <- paste0("NoWol.", fdr.col)
      fc.col <- paste0("NoWol.", fc.col)
    }
    signif.rows <- ss.model.results[,fdr.col] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    signif.rows[failed.rows] <- F
    if(sum(failed.rows) > 0) {
      cat(sum(failed.rows),"features failed to converge when testing Treatment fixed effect ANOVA\n")
      cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
    }
    cat(sum(signif.rows),"features have significant Treatment effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
    cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n\n")
    cat("Distribution of Treatment Effects:\n")
    cat("Range:", range(ss.model.results[,fc.col]), "\n")
    cat("Box Stats:", boxplot.stats(ss.model.results[,fc.col])$stats,"\n")
    cat("\n")
  }
  
  # Store this table in ss.results list
  ss.results[[ss]] <- ss.model.results
  
  # Output line means ONLY for non-permutation analysis
  if(is.na(my.args$PERM)) {
    # UNSCALE HERE
    for(j in 1:ncol(ss.line.means)) {
      ss.line.means[,j] <- (ss.line.means[,j] * ss.sd) + ss.center
    }
    
    # Create output table, first column is gene ID, each subsequent column is group mean expr fitted from regression
    line.mean.table <- cbind(GENE=row.names(ss.line.means), ss.line.means)
    
    # Remove genes that were not significantly heritable here
    test.cols <- paste0(c("Line", "Sex.Line", "Trmt.Line", "Sex.Trmt.Line"), ".PVal.FDR")
    if(my.args$WOLADJ) {
      test.cols <- paste0("WolAdj.", test.cols)
    } else {
      test.cols <- paste0("NoWol.", test.cols)
    }
    test.cols <- intersect(test.cols, colnames(ss.model.results))
    stopifnot(length(test.cols) > 0)
    
    if(length(test.cols) == 1) {
      line.mean.rows <- row.names(ss.model.results)[ss.model.results[,test.cols] <= my.args$FDR]
    } else {
      line.mean.rows <- row.names(ss.model.results)[apply(ss.model.results[,test.cols], 1, function(q){any(q <= my.args$FDR)})]
    }
    cat("Outputting line means for", length(line.mean.rows), "features with FDR <", my.args$FDR, "in any of these columns:\n")
    cat(test.cols, "\n")
    
    line.mean.table <- line.mean.table[line.mean.rows,]
  
    # Output the table to disk, overwrite if it exists
    line.mean.file <- paste0(my.args$OUTSTUB, if(my.args$WOLADJ){"_WolAdj_"}else{"_NoWol_"}, ss, "_line_means.txt")
    cat("Writing line means for", ss, "samples to:", line.mean.file, "\n")
    if(file.exists(line.mean.file)) {
      cat("WARNING: Overwriting existing file\n")
    }
    write.table(line.mean.table, line.mean.file, sep="\t", row.names=F, quote=F)
  }
  
  cat("\nFinished linear mixed model analysis for", ss, "samples.\n\n\n")
}

# First put subset ID in front of each column name
# Also fill in values for skipped transcripts
# TO DO: This should probably just be incorporated in the loop above
all.features <- row.names(expr.table)
for(ss in names(ss.results)) {
  missing.features <- setdiff(all.features, row.names(ss.results[[ss]]))
  if(length(missing.features) > 0) {
    # Set all PercVar, Log2FC, and H2 columns to 0
    zero.cols <- c(grep("Perc[.]Var$", colnames(ss.results[[ss]])), grep("Log2FC$", colnames(ss.results[[ss]])), grep("H2$", colnames(ss.results[[ss]])))
    ss.results[[ss]][missing.features,zero.cols] <- 0
    # Set all PVal/FDR columns to 1
    pval.cols <- grep("PVal", colnames(ss.results[[ss]]))
    ss.results[[ss]][missing.features,pval.cols] <- 1
  }
  colnames(ss.results[[ss]]) <- paste(ss, colnames(ss.results[[ss]]), sep=".")
}

# Combine all subset result tables
all.ss.results <- data.frame(row.names=all.features)
for(ss in names(ss.results)) {
  all.ss.results <- cbind(all.ss.results, ss.results[[ss]][all.features,])
}

# Fast reload to skip re-running this step during testing:
# all.ss.results <- read.table(my.args$H2FILE, header=T, sep="\t", row.names=1)

# Summarize overall, depending on model type:
cat("\nSummary of results:\n\n")

# Built up the list of appropriate FDR columns to use
fdr.cols <- c(LINE="Line.PVal.FDR")
if(my.args$POOLED) {
  fdr.cols <- c(fdr.cols, 'SEX:LINE'="Sex.Line.PVal.FDR")
}
if(use.trmt) {
  fdr.cols <- c(fdr.cols, 'TRMT:LINE'="Trmt.Line.PVal.FDR")
}
if(my.args$POOLED & use.trmt) {
  fdr.cols <- c(fdr.cols, 'SEX:TRMT:LINE'="Sex.Trmt.Line.PVal.FDR")
}
fdr.col.names <- names(fdr.cols)

if(my.args$WOLADJ) {
  fdr.cols <- paste0("WolAdj.", fdr.cols)
} else {
  fdr.cols <- paste0("NoWol.", fdr.cols)
}
names(fdr.cols) <- fdr.col.names

if(length(ss.results) == 1) {
  # Summarization of pooled or single-sex model
  fdr.cols <- paste(names(ss.results), fdr.cols, sep=".")
  names(fdr.cols) <- fdr.col.names
  
  for(term in names(fdr.cols)) {
    signif.term <- all.ss.results[,fdr.cols[term]] <= my.args$FDR
    cat(sum(signif.term), "features with significant", term, "effects =", round(mean(signif.term)*100), "%\n")
  }
  
  if(length(fdr.cols) > 1) {
    signif.all <- apply(all.ss.results[,fdr.cols], 1, function(x){all(x <= my.args$FDR)})
    signif.any <- apply(all.ss.results[,fdr.cols], 1, function(x){any(x <= my.args$FDR)})
    cat(sum(signif.all), "features are significant for ALL terms =", round(mean(signif.all)*100), "%\n")
    cat(sum(signif.any), "features are significant for at least ONE term =", round(mean(signif.any)*100), "%\n")
  }
  cat("\n")
  
  # Split by known/novel genes if appropriate:
  if(split.known.novel) {
    cat("Breakdown for KNOWN genes:\n")
    for(term in names(fdr.cols)) {
      signif.term <- all.ss.results[known.rows,fdr.cols[term]] <= my.args$FDR
      cat(sum(signif.term), "known genes with significant", term, "effects =", round(mean(signif.term)*100), "%\n")
    }
    if(length(fdr.cols) > 1) {
      cat(sum(signif.all[known.rows]), "known genes are significant for ALL terms =", round(mean(signif.all[known.rows])*100), "%\n")
      cat(sum(signif.any[known.rows]), "known genes are significant for at least ONE term =", round(mean(signif.any[known.rows])*100), "%\n")
    }
    cat("\n")
    
    cat("Breakdown for NOVEL genes:\n")
    for(term in names(fdr.cols)) {
      signif.term <- all.ss.results[novel.rows,fdr.cols[term]] <= my.args$FDR
      cat(sum(signif.term), "novel genes with significant", term, "effects =", round(mean(signif.term)*100), "%\n")
    }
    if(length(fdr.cols) > 1) {
      cat(sum(signif.all[novel.rows]), "novel genes are significant for ALL terms =", round(mean(signif.all[novel.rows])*100), "%\n")
      cat(sum(signif.any[novel.rows]), "novel genes are significant for at least ONE term =", round(mean(signif.any[novel.rows])*100), "%\n")
    }
    cat("\n")    
  }
} else {
  # Summarization of within-sex models for both sexes
  for(term in names(fdr.cols)) {
    term.fdr.cols <- paste(names(ss.results), fdr.cols[term], sep=".")
    signif.any.sex <- apply(all.ss.results[,term.fdr.cols], 1, function(x){any(x <= my.args$FDR)})
    signif.all.sex <- apply(all.ss.results[,term.fdr.cols], 1, function(x){all(x <= my.args$FDR)})
    cat(sum(signif.any.sex), "features have significant", term, "effects in at least one sex =", round(mean(signif.any.sex)*100), "%\n")
    cat(sum(signif.all.sex), "features have significant", term, "effects in BOTH sexes =", round(mean(signif.all.sex)*100), "%\n")
    cat("\n")
  }
  
  if(length(fdr.cols) > 1) {
    all.fdr.cols <- paste(rep(names(ss.results), each=length(fdr.cols)), rep(fdr.cols, times=length(ss.results)), sep=".")
    signif.any <- apply(all.ss.results[,all.fdr.cols], 1, function(x){any(x <= my.args$FDR)})
    cat(sum(signif.any), "features are significant for at least one term in at least one sex =", round(mean(signif.any)*100), "%\n")
    cat("\n")
  }
  
  # Split by known/novel genes if appropriate:
  if(split.known.novel) {
    cat("Breakdown for KNOWN genes:\n")
    for(term in names(fdr.cols)) {
      term.fdr.cols <- paste(names(ss.results), fdr.cols[term], sep=".")
      signif.any.sex <- apply(all.ss.results[known.rows,term.fdr.cols], 1, function(x){any(x <= my.args$FDR)})
      signif.all.sex <- apply(all.ss.results[known.rows,term.fdr.cols], 1, function(x){all(x <= my.args$FDR)})
      cat(sum(signif.any.sex), "known genes have significant", term, "effects in at least one sex =", round(mean(signif.any.sex)*100), "%\n")
      cat(sum(signif.all.sex), "known genes have significant", term, "effects in BOTH sexes =", round(mean(signif.all.sex)*100), "%\n")
      cat("\n")
    }
    
    if(length(fdr.cols) > 1) {
      cat(sum(signif.any[known.rows]), "known genes are significant for at least one term in at least one sex =", round(mean(signif.any[known.rows])*100), "%\n")
      cat("\n")
    }
    
    cat("Breakdown for NOVEL genes:\n")
    for(term in names(fdr.cols)) {
      term.fdr.cols <- paste(names(ss.results), fdr.cols[term], sep=".")
      signif.any.sex <- apply(all.ss.results[novel.rows,term.fdr.cols], 1, function(x){any(x <= my.args$FDR)})
      signif.all.sex <- apply(all.ss.results[novel.rows,term.fdr.cols], 1, function(x){all(x <= my.args$FDR)})
      cat(sum(signif.any.sex), "novel genes have significant", term, "effects in at least one sex =", round(mean(signif.any.sex)*100), "%\n")
      cat(sum(signif.all.sex), "novel genes have significant", term, "effects in BOTH sexes =", round(mean(signif.all.sex)*100), "%\n")
      cat("\n")
    }
    
    if(length(fdr.cols) > 1) {
      cat(sum(signif.any[novel.rows]), "novel genes are significant for at least one term in at least one sex =", round(mean(signif.any[novel.rows])*100), "%\n")
      cat("\n")
    }
  }
}

# Output full result table to a file
cat("\nWriting all model results to:", my.args$H2FILE, "\n")
write.table(cbind(GENE=row.names(all.ss.results), all.ss.results), my.args$H2FILE, sep="\t", row.names=F, quote=F)

cat("\n\nScript completed successfully!\n\n")
print(proc.time())
quit("no")
