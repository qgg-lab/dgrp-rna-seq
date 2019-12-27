#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 10/6/16
# MAJOR UPDATE - 2/25/18
# Corrected issues with line mean/diff calculations
#
# gen_var_model.R
# (Adapted from/Replaces assess_heritability.R)
#
# Goal: Load in normalized data and apply linear mixed model to assess H^2 for each gene/feature
# This version of the script defines the model in a more flexible way
# 
# Usage:
# Rscript assess_heritability.R [OPTIONS] EXPR=subdir/data_[rpm|fpkm|rle].txt
#  EXPR=  The file containing expression values to assess heritability on
#  SAMPLES= The table with sample information (defaults to sample_master_table.txt)
#  LINEINFO=  The file containing line info (Wolbachia, Inversions, PCs), 
#          defaults to ~ljeveret/Resources/DGRP/Freeze2/line_adjust_info.txt
#  FILTER=  Flag(s) to filter genes on. Multiple flags should be comma-seperated, no spaces
#           Default is to leave it blank, no filtering on flag column (will still get passed through to output tables)
#           FILTER=LOW filters out just the genes with "LOW" in FLAG column
#           FILTER=LOW,RARE filters out genes with "LOW" or "RARE" in FLAG column
#  WOLADJ=  TRUE/FALSE (Default=TRUE), adjust for Wolbachia effects before computing H^2 and line means
#  EQTL=   TRUE/FALSE (Default=FALSE), Incorporates Major Inversions and genetic PC correction factors in the model as well
#     If TRUE, these will also be removed from line means
#  POOLED=  TRUE/FALSE (Default=FALSE), Run a pooled model instead of within-sex models
#  TRMT=  Name of an additional column in SAMPLES to use as Treatment variable in model (default=None)
#  CORES= The number of CPUs to use for parallelization (defaults to 1)
#   NOTE: If submitting your job through SLURM, make sure to set -c option as well
#  NQ=    TRUE/FALSE (Default=FALSE), apply Normal Quantile transformation before running H^2
#  LINEREG=  File with gene x line scores to regress over BEFORE fitting main genetic variance model. 
#           If blank, no adjustment for this is performed (Default)
#           But if a file is specified, then modeling will be performed ONLY on the intersect set of genes
#  LRLOG= Convert all LINEREG values to log2 scale before fitting models (defaults to T)
#  TAG=   Extra tag added to output file names, 
#         default is leave blank unless LINEREG is set, in which case add "LR"
#  PERM=  If specified, permute the sample labels within each group (sex)
#         This should be provided as an integer value, which is used to both seed the RNG
#         And give the output a unique ID (output will have _PermN_ designation code)
#         All output goes in Perm/ subdir of Out Dir as well
#  OUTDIR= Directory to output, default is same directory where EXPR file is.
#  MEANS= T/F, Output line means? Defaults to TRUE unless PERM is set
#  FDR=   FDR cut-off (for reporting counts and outputtin line means), default = 0.05
#  PSEUDO=  Pseudo-count used to convert RPM or FPKM values to log scale, default = 0.001
#  LRPSEUDO=  Pseudo-count used specifically to convert the LINEREG values if doing Log2 conversion
#             Defaults to same as PSEUDO, which should be appropriate for rpm type values (DNA Adjustment), but NOT variant rate
#
# Automatically puts all output in the same subdir as the data file, with similar basename
# Normally takes log2 of the data (with pseudo count), UNLESS data file ends in *_rle
# Automatically uses sample_master_table.txt for sample information
#
# TO DO: Should Pooled Treatment model also include SEX:TRMT term?
#
# TO DO: Make option to drop any lines that are "incomplete" (< 2 reps within sex, pooled version already does this)
#
# TO DO: Make interface more flexible:
#   Let each param have optional leading dashes, and capitalization doesn't matter
#


# --- INITIALIZATION --- #

options(stringsAsFactors=F)

# NOTE: With Levene's test no longer done in this script, this dependency is no longer needed!
# Standard package needed for Levene's test
# If missing, run the following in same R version:
# install.packages("car")
# library(car)

# Main package for the linear mixed model
# If not installed, run:
# install.packages(c("lme4","lmerTest","doMC"))
# NOTE: Only need to load lmerTest - it loads what it needs from lme4 or has the same functionality
# library(lme4)
library(lmerTest)
# This one is needed for parallelization later on
library(doMC)


# -- Fixed Parameters -- #

# Terms used in eQTL mode
eqtl.inv <- c("In2Lt","In2RNS","In3RP","In3RK","In3RMo")
eqtl.pc <- paste0("PC",1:10)


# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript assess_heritability.R [OPTIONS] EXPR=subdir/data_[rpm|fpkm|rle].txt"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING (Baseline):
# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")
# my.args <- c("EXPR=microbiome/combined_microbe_filtered_rpm.txt","OUTDIR=microbiome/genVar/","LINEINFO=~/Resources/DGRP/Freeze2/line_adjust_info.txt")
# my.args <- c("POOLED=TRUE","EXPR=microbiome/combined_microbe_filtered_rpm.txt","OUTDIR=microbiome/genVar/","EQTL=TRUE","LINEINFO=~/Resources/DGRP/Freeze2/line_adjust_info.txt")
# INTERACTIVE TESTING (3WK):
# setwd("~/Projects/DGRP_3WK_RNAseq_Post/")
# my.args <- c("POOLED=TRUE","WOLADJ=FALSE","TRMT=AGE","EXPR=microbiome/combined_microbe_special_rpm.txt","OUTDIR=$MYPATH/genVar/","LINEINFO=~/Resources/DGRP/Freeze2/line_adjust_info.txt")
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

setDefault("LINEINFO","~ljeveret/Resources/DGRP/Freeze2/line_adjust_info.txt")

setDefault("FILTER","")

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

setDefault("EQTL",F)
# If EQTL is "0" or "1" map these to "FALSE" and "TRUE"
if(is.character(my.args$EQTL) & (my.args$EQTL == "0")){my.args$EQTL <- F}
if(is.character(my.args$EQTL) & (my.args$EQTL == "1")){my.args$EQTL <- T}
# Convert any other type to logical
if(!is.logical(my.args$EQTL)){my.args$EQTL <- as.logical(my.args$EQTL)}
# If conversion produced NA, just set to T
if(is.na(my.args$EQTL)) {
  cat("WARNING: Did not recognize value of EQTL parameter, setting to F by default.\n")
  my.args$WOLADJ <- F
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

setDefault("LINEREG","")

setDefault("LRLOG",T)
# If LRLOG is "0" or "1" map these to "FALSE" and "TRUE"
if(is.character(my.args$LRLOG)) {
  if(my.args$LRLOG == "0"){my.args$LRLOG <- F}
  if(my.args$LRLOG == "1"){my.args$LRLOG <- T}
}
# Convert any other type to logical
if(!is.logical(my.args$LRLOG)){my.args$LRLOG <- as.logical(my.args$LRLOG)}
# If conversion produced NA, just set to T
if(is.na(my.args$LRLOG)) {
  cat("WARNING: Did not recognize value of LRLOG parameter, setting to T by default.\n")
  my.args$LRLOG <- T
}

setDefault("TAG", if(my.args$LINEREG == ""){""}else{"LR"})

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

setDefault("MEANS",is.na(my.args$PERM))

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

setDefault("PSEUDO", 0.001)
if(!is.numeric(my.args$PSEUDO)) {
  my.args$PSEUDO <- as.numeric(my.args$PSEUDO)
}

setDefault("LRPSEUDO", my.args$PSEUDO)
if(!is.numeric(my.args$LRPSEUDO)) {
  my.args$LRPSEUDO <- as.numeric(my.args$LRPSEUDO)
}

setDefault("FDR", 0.05)

# OUTDIR default: Extract path to EXPR input file, use it as the output dir, append Perm if doing permutation analysis
setDefault("OUTDIR", paste0(dirname(my.args$EXPR), "/", if(!is.na(my.args$PERM)){"Perm/"}))
if(!grepl("[/]$", my.args$OUTDIR)) {
  my.args$OUTDIR <- paste0(my.args$OUTDIR, "/")
}
# Make sure this subdir exists
if(!file.exists(my.args$OUTDIR)) {
  dir.create(path=my.args$OUTDIR, recursive=T, showWarnings = F)
}

# Base name of COUNTS file comes from shaving off subdirectory part of path, and trailing ".txt"
my.args$BASENAME <- sub("[.]txt$", "", basename(my.args$EXPR))

# OUTSTUB will be OUTDIR/BASENAME(_OptionFlags)
# This is used to determine output file paths
my.args$OUTSTUB <- paste0(my.args$OUTDIR, my.args$BASENAME)
if(my.args$NQ) {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_NQ")
}
if(my.args$TAG != "") {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_", my.args$TAG)
}
if(my.args$EQTL) {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_eQTL")
} else if(my.args$WOLADJ) {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_WolAdj")
}
if(!my.args$WOLADJ) {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_NoWol")
}
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
if(grepl("[_.]counts$", my.args$BASENAME)) {
  my.args$TYPE <- "counts"
  # my.args$BASENAME <- sub("[_.]counts$", "", my.args$BASENAME)
  cat("WARNING: Input file appears to be raw count data, will transform to log2 scale, but this script is intended for normalized expression (rpm, fpkm, or rle)!\n")
} else if(grepl("[_.]rpm$", my.args$BASENAME)) {
  my.args$TYPE <- "rpm"
  # my.args$BASENAME <- sub("[_.]rpm$", "", my.args$BASENAME)
  cat("Treating input expression file as reads per million (rpm), will transform to log2 scale before fitting lme.\n")
} else if(grepl("[_.]fpkm$", my.args$BASENAME)) {
  my.args$TYPE <- "fpkm"
  # my.args$BASENAME <- sub("[_.]fpkm$", "", my.args$BASENAME)
  cat("Treating input expression file as features per kb per million (fpkm), will transform to log2 scale before fitting lme.\n")
} else if(grepl("[_.]rle$", my.args$BASENAME)) {
  my.args$TYPE <- "rle"
  # my.args$BASENAME <- sub("[_.]rle$", "", my.args$BASENAME)
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
# If first column == FLAG, use that to filter the table (optional), then save in a separate named vector
# NOTE: This no longer filters OTHER flags (MT,RIBO,etc.) automatically
#  - those should be removed by normalization script, but if not use the FILTER option to remove them
dropped.features <- 0
if(colnames(expr.table)[1] == "FLAG") {
  if(my.args$FILTER!="") {
    my.args$FILTER <- unlist(strsplit(my.args$FILTER, split=","))
    drop.rows <- expr.table$FLAG %in% my.args$FILTER
    dropped.features <- sum(drop.rows)
    expr.table <- expr.table[!drop.rows,]
  }
  gene.flags <- expr.table$FLAG
  names(gene.flags) <- row.names(expr.table)
  expr.table <- expr.table[,2:ncol(expr.table)]
} else {
  # Fill in flag column with OK for all
  gene.flags <- rep("OK", times=nrow(expr.table))
  names(gene.flags) <- row.names(expr.table)
}
# Remove 'X' that R automatically puts in front of column names
colnames(expr.table) <- sub("^X", "", colnames(expr.table))
features <- row.names(expr.table)
cat("Loaded counts for",length(features),"features across",ncol(expr.table),"samples from",my.args$EXPR,"\n")
if(dropped.features > 0) {
  cat("(Dropped",dropped.features,"features based on FLAG =",paste(my.args$FILTER, collapse=" or "),")\n")
}
cat("Breakdown of flag counts (after filtering):\n")
table(gene.flags)

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
# Remove DGRP_ from row names
row.names(line.table) <- sub("^DGRP[_]","",row.names(line.table))
# Make sure all lines in sample.table are in line.table
stopifnot(all(as.character(sample.table$LINE) %in% row.names(line.table)))
# Subset to just the lines in sample.table
line.table <- line.table[unique(as.character(sample.table$LINE)),]
cat("Loaded info on",nrow(line.table),"lines from",my.args$LINEINFO,"\n")
cat(sum(line.table$Wolbachia),"lines have Wolbachia.\n")
sample.table$WOL <- line.table[as.character(sample.table$LINE),"Wolbachia"]
cat(round(mean(sample.table$WOL)*100),"% of samples have Wolbachia.\n")
cat("\n")

# TO DO: Script needs to be smarter about cases where WOL or any given INV only has a single group
# Need to check for these cases and put the relevant factors in a "drop.terms" vector, drop below

# For eQTL mode, must also contain Inversions and 10 PCs - copy these over to sample.table as well
# NOTE: In Wen's script he converts Inversions to factors, but also may have started with them in 0/1/2 format
# TO DO: Talk to Wen about whether each Inv should be modeled as having a linear effect (INV effect 2x that of HET effect?)
if(my.args$EQTL) {
  inv.cols <- grep("^I[nN]", colnames(line.table))
  if(length(inv.cols) == 0) {
    stop(my.args$LINEINFO, " did not contain any Inversion columns, cannot run in EQTL mode.")
  }
  # Remove [_] from Inv cols
  colnames(line.table)[inv.cols] <- gsub("[_]","",colnames(line.table)[inv.cols])
  # Make sure the req'd Inversions are present:
  if(!all(eqtl.inv %in% colnames(line.table))) {
    stop(my.args$LINEINFO, " is missing one or more Inv columns for EQTL correction:", paste(setdiff(eqtl.inv, colnames(line.table)), collapse=", "))
  }
  # Make sure each Inversion column has only 3 different values, convert to factor
  inv.groups <- c("ST","HET","INV")
  for(j in eqtl.inv) {
    if(!all(line.table[,j] %in% inv.groups)) {
      stop(my.args$LINEINFO, " contains non-standard values for ", j, "column")
    }
    line.table[,j] <- factor(line.table[,j], levels=intersect(inv.groups, unique(line.table[,j])))
  }
  # Make sure the req'd PCs are present:
  eqtl.pc <- paste0("PC", 1:10)
  if(!all(eqtl.pc %in% colnames(line.table))) {
    stop(my.args$LINEINFO, " is missing one or more PC columns for EQTL correction:", paste(setdiff(eqtl.pc, colnames(line.table)), collapse=", "))
  }
  
  sample.table <- cbind(sample.table, line.table[as.character(sample.table$LINE),c(eqtl.inv,eqtl.pc)])
}

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
  # THIS STEP IS DEPRECATED
  # TO DO: ADD OPTION TO AUTOMATICALLY DROP LOW AND/OR RARE FEATURES HERE (OR ABOVE?) BASED ON FLAGS
  # SHOULD FIRST CONFIRM THAT MOST OF THESE FEATURES DO NOT GIVE SIGNIF GENVAR
  #---------
  # min.rep.expr.perc <- apply(expr.table, 1, function(x){mean(x > 0)})
  #
  # non.rep.features <- names(min.rep.expr.perc)[min.rep.expr.perc < 0.05]
  #
  # cat("Removing", length(non.rep.features), "features with expression in <",my.args$FDR*100,"% of all samples.\n")
  # keep.features <- setdiff(row.names(expr.table), non.rep.features)
  # expr.table <- expr.table[keep.features,]
  
  cat("Converting to log2 scale, with pseudo-count =", my.args$PSEUDO, "\n")
  expr.table <- log2(expr.table + my.args$PSEUDO)
} else {
  cat("No filtering or log2 transformation applied, data is assumed to be log2 scale already (rle type).\n")
}

# If specified, load the LINEREG file
linereg.table <- NULL
if(my.args$LINEREG != "") {
  stopifnot(file.exists(my.args$LINEREG))
  cat("Loading line regression data from:", my.args$LINEREG, "\n")
  linereg.table <- read.table(my.args$LINEREG, header=T, sep="\t", row.names=1)
  cat("Loaded line regression data for", nrow(linereg.table), "features in", length(setdiff(colnames(linereg.table), "TOTAL")), "lines.\n")
  keep.genes <- row.names(expr.table) %in% row.names(linereg.table)
  if(!all(keep.genes)) {
    cat("WARNING: Missing line regression data for", sum(!keep.genes), "features.\n")
    cat("These features will be dropped and H2 model will only run for", sum(keep.genes), "features.\n")
    expr.table <- expr.table[keep.genes,]
  }
  # Drop line regression data for genes missing from expr table (don't need to warn on this one)
  stopifnot(all(row.names(expr.table) %in% row.names(linereg.table)))
  linereg.table <- linereg.table[row.names(expr.table),]
  # Replace leading "X" with "line_" on column names
  colnames(linereg.table) <- sub("^X", "line_", colnames(linereg.table))
  # Remove FLAG column if present
  if("FLAG" %in% colnames(linereg.table)) {
    linereg.table <- linereg.table[,setdiff(colnames(linereg.table), "FLAG")]
  }
  # Convert to log2 scale (optional)
  if(my.args$LRLOG) {
    cat("Converting line regression data to log2 scale with pseudo-count =", my.args$LRPSEUDO, "\n")
    linereg.table <- log2(linereg.table + my.args$LRPSEUDO)
  }
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


# --- Core Sub-Routines --- #

# NOTE: Splitting out the individual model-fitting steps into sub-routines seemed to 
# cause problems because passing a complete model result back from a function appears to lose some critical data
# This is a flaw in the lmer module design as best I can tell

# Functions to extract just the fixed or random effect terms, or just the primary or interaction terms from a vector of model terms
getFixedTerms <- function(my.terms) {
  grep("^[(]1[|].*[)]$", my.terms, invert=T, value=T)
}

getRandomTerms <- function(my.terms) {
  grep("^[(]1[|].*[)]$", my.terms, value=T)
}

getPrimaryTerms <- function(my.terms) {
  grep("[:]", my.terms, invert=T, value=T)
}

getIxnTerms <- function(my.terms) {
  grep("[:]", my.terms, value=T)
}

# Function to extract the names of all columns that need to be present in the input based on the vector of model terms
getInputTerms <- function(my.terms) {
  my.terms.split <- unlist(strsplit(my.terms, split="[()|:]"))
  return(setdiff(unique(my.terms.split), c("1","")))
}

# Function to build complete model string based on model terms (predicted value is "EXPR" by default)
getModelString <- function(my.terms, prediction="EXPR") {
  paste0(prediction, " ~ ", paste(my.terms, collapse=" + "))
}

# !!! IMPT !!!
# TO DO: This should convert all model term columns to factors!
# This function combines one expression vector with all model data
getModelInput <- function(
  i,  # Gene index to extract a model input table for
  my.expr,  # Expression matrix with columns = genes, samples = rows (should also be scaled and normalized already)
  my.samples,  # data frame with samples = rows, and all columns indicated by model.inputs and group.col
  model.inputs=c("SEX","WOL","LINE","TRMT"),  # Model input columns to include in return table
  group.col="GROUP",  # Column with group names used for line means table
  line.reg=NULL  # Specifiy the line reg table here (same dimensionality as my.expr) when you want to include LINEREG as a model correction factor
) {
  stopifnot(nrow(my.expr)==nrow(my.samples))
  stopifnot(row.names(my.expr)==row.names(my.samples))
  my.model.input <- cbind(EXPR=my.expr[,i], my.samples[,c(model.inputs,group.col)])
  # Convert the model inputs to factors (except for PC columns!)
  convert.cols <- setdiff(model.inputs, grep("^PC[0-9]",model.inputs,value=T))
  # Also skipping Inversion columns here b/c they should already be converted to factors (and we may want to try fitting them as additive)
  convert.cols <- setdiff(convert.cols, grep("^I[nN][2-4XY]",convert.cols,value=T))
  for(j in convert.cols) {
    my.model.input[,j] <- factor(my.model.input[,j])
  }
  if(!is.null(line.reg)) {
    # Extract variant rate data if provided
    stopifnot(all(dim(my.expr)==dim(line.reg)))
    stopifnot(all(row.names(my.expr)==row.names(line.reg)))
    stopifnot(all(colnames(my.expr)==colnames(line.reg)))
    my.model.input[,"LINEREG"] <- line.reg[,i]
  }
  return(my.model.input)
}

# This function takes a single model input data frame and performs ALL modeling steps
# my.model.input should be a data frame as returned by getModelInput above
# my.model.terms should be a character vector with all terms used in the core model
# Set verbose=T to see progress reporting - DO NOT USE FOR FULL/PARALLEL PROCESSING!
# All other params currently pulled from global my.args
# Return type is a list with two members:
# $model.results = vector of all parameters for model result table (except adjusted p-values)
# $line.means = vector of line means
# NOTE: This function does not handle SCALING any model params or line means back to the original distributions
#       The model.input is taken at face value, without any scaling data
#       This function also does not handle p-value adjustments, since it only sees one gene expression vector at a time
# TO DO: Could set all parameters explicitly, rather than relying on globals
#     (only a few are actually used here: my.args$POOLED, my.args$WOLADJ, use.trmt)
modelGeneticVar <- function(my.model.input, my.model.terms, verbose=F) {
  # This vector will store named model result values:
  my.results <- c()
  
  # If there's a LINEREG column present, regress out these effects first
  if("LINEREG" %in% colnames(my.model.input)) {
    # Regress out by fitting a standard LM first
    # Residuals will be the new expr for main model
    # NOTE: We *COULD* try incorporating these directly in the main model, but I think this is the better way to do things
    # This pulls out the maximum amount of variance attributable to alignment bias BEFORE fitting the main model of interest
    if(verbose) cat("Correcting for gene level variation rate.\n")
    # TO DO: For POOLED or TRMT model, could include ixn terms here?
    #   This would be for genes that are only expressed in one sex or one treatment group
    #   But I don't think this is necessary - we don't really care about this correction in either of those contexts yet
    #   Would also need to add back the primary SEX and/or TRMT effects here
    linereg.lm <- lm(EXPR ~ LINEREG, data=my.model.input)
    
    # Compute Coefficient and Avg Effect Size (This can be scaled back later)
    my.results["LINEREG.Coeff"] <- linereg.lm$coefficients["LINEREG"]
    my.results["LINEREG.Log2FC"] <- mean(linereg.lm$coefficients["LINEREG"] * my.model.input$LINEREG)
    # These can be NA if this gene has no recorded variants, so LINEREG is 0 for all
    if(is.na(my.results["LINEREG.Coeff"])) my.results["LINEREG.Coeff"] <- 0
    if(is.na(my.results["LINEREG.Log2FC"])) my.results["LINEREG.Log2FC"] <- 0
    
    # Compute the corrected expression values
    # Add intercept back into each column of resid to keep similar scale
    adj.expr <- resid(linereg.lm) + linereg.lm$coefficients["(Intercept)"]
    
    orig.var <- var(my.model.input$EXPR)
    adj.var <- var(adj.expr)
    my.results["LINEREG.Perc.Var"] <- (orig.var - adj.var) / orig.var
    
    # Use ANOVA to get P-value
    my.results["LINEREG.PVal"] <- anova(linereg.lm)$Pr[1]
    # If NA, set to P-value=1
    if(is.na(my.results["LINEREG.PVal"])) my.results["LINEREG.PVal"] <- 1
    
    # Replace the original expression data with the adjusted expression data
    my.model.input$EXPR <- adj.expr
  }
  
  # Fit the flexible model
  my.formula <- formula(getModelString(my.model.terms))
  if(verbose) {
    cat("Fitting the following model:\n")
    print(my.formula)
    cat("\n")
  }
  my.lmer <- lmer(my.formula, data=my.model.input)
  
  # As long as there is 1 or more fixed effect terms, perform anova tests, compute Log2FC, and compute Perc.Var
  fTerms <- getFixedTerms(my.model.terms)
  if(length(fTerms) > 0) {
    
    # Extract P-values from ANOVA test
    my.anova <- anova(my.lmer)
    stopifnot(all(fTerms %in% row.names(my.anova)))
    if("Pr(>F)" %in% colnames(my.anova)) {
      fixed.pvals <- my.anova$Pr
    } else {
      fixed.pvals <- rep(as.numeric(NA), times=nrow(my.anova))
    }
    names(fixed.pvals) <- row.names(my.anova)
    fixed.pvals <- fixed.pvals[fTerms]
    names(fixed.pvals) <- paste0(names(fixed.pvals), ".PVal")
    
    # Identify fixed effect terms with multiple (>2) levels
    # This list structure matches each fixed effect terms to the complete set of model fixed effects
    fTerm.effects <- lapply(getPrimaryTerms(fTerms), function(ft){
      if(class(my.model.input[,ft])=="factor") {
        ft.levels <- levels(my.model.input[,ft])
        stopifnot(length(ft.levels) > 1)
        ft.cols <- paste0(ft,ft.levels[-1])
        if(length(ft.cols) > 1) {
          names(ft.cols) <- paste(ft,ft.levels[-1],sep=".")
        } else {
          names(ft.cols) <- ft
        }
      } else {
        ft.cols <- ft
        names(ft.cols) <- ft
      }
      return(ft.cols)
    })
    names(fTerm.effects) <- getPrimaryTerms(fTerms)
    for(ft in getIxnTerms(fTerms)) {
      ft.primaries <- unlist(strsplit(ft, ":", fixed=T))
      ft.prim.levels <- lapply(ft.primaries, function(x){fTerm.effects[[x]]})
      ft.ixn.levels <- apply(expand.grid(ft.prim.levels), 1, paste, collapse=":")
      ft.pl.names <- lapply(ft.primaries, function(x){names(fTerm.effects[[x]])})
      names(ft.ixn.levels) <- apply(expand.grid(ft.pl.names), 1, paste, collapse=":")
      fTerm.effects[[ft]] <- ft.ixn.levels
    }
    
    # Extract the Log2FC values - keep in list structure for now
    my.fixef <- fixef(my.lmer)
    fixed.log2fc <- lapply(fTerm.effects, function(ft.eff){
      ft.log2fc <- my.fixef[ft.eff]
      names(ft.log2fc) <- paste0(names(ft.eff), ".Log2FC")
      return(ft.log2fc)
    })
    
    # Compute % Var here by taking the variance of actual expression vs fixed-effect corrected expr var
    # Total variance of uncorrected expression values
    orig.var <- var(my.model.input$EXPR)
    # First compute the matrix of fixed effects by all samples
    my.model <- model.matrix(my.lmer)
    my.model <- t(apply(my.model, 1, function(x){x * my.fixef}))
    stopifnot(all(row.names(my.model) == row.names(my.model.input)))
    # Collapse my.model columns (1 per term level) into columns for each term (add together columns for levels of same term)
    effect.matrix <- matrix(as.numeric(NA), ncol=length(fTerms), nrow=nrow(my.model), dimnames=list(row.names(my.model),fTerms))
    for(ft in fTerms) {
      model.cols <- fTerm.effects[[ft]]
      stopifnot(length(model.cols) > 0)
      if(length(model.cols) == 1) {
        effect.matrix[,ft] <- my.model[,model.cols]
      } else {
        effect.matrix[,ft] <- apply(my.model[,model.cols], 1, sum)
      }
    }
    stopifnot(sum(is.na(effect.matrix))==0)
    
    fixed.perc.var <- c()
    for(j in colnames(effect.matrix)) {
      if(grepl("[:]", j)) {
        # TO DO: Confirm this still works on new version?
        # For fixed effect interaction terms, compare correction for the individual fixed effects vs the ixn
        j.components <- unlist(strsplit(j, "[:]"))
        j.indiv.only <- my.model.input$EXPR
        for(jc in j.components) {
          j.indiv.only <- j.indiv.only - effect.matrix[,jc]
        }
        j.ixn <- j.indiv.only - effect.matrix[,j]
        j.indiv.var <- var(j.indiv.only)
        j.ixn.var <- var(j.ixn)
        fixed.perc.var[j] <- (j.indiv.var - j.ixn.var) / j.indiv.var
      } else {
        # For primary fixed effect, just compare original expression variance to variance after correcting with this term
        j.sub.var <- var(my.model.input$EXPR - effect.matrix[,j])
        fixed.perc.var[j] <- (orig.var - j.sub.var) / orig.var
      }
    }
    names(fixed.perc.var) <- paste0(names(fixed.perc.var), ".Perc.Var")
    
    # Linearize so that results have the form Term1.PVal, Term1.(Levels).Log2FC, Term1.PercVar, Term2.PVal, ...
    for(ft in fTerms) {
      my.results <- c(my.results, fixed.pvals[paste0(ft,".PVal")], fixed.log2fc[[ft]], fixed.perc.var[paste0(ft,".Perc.Var")])
    }
  }
  
  # As long as there is 1 or more random effect term, perform rand tests, compute H2
  if(length(getRandomTerms(my.model.terms)) > 0) {
    # Compute varcorr of all terms
    my.varcorr <- as.data.frame(VarCorr(my.lmer))
    my.varcorr.names <- my.varcorr$grp
    my.varcorr <- my.varcorr$vcov
    names(my.varcorr) <- my.varcorr.names
    
    # H2 is every term involving LINE vs sum of all variance
    line.terms <- grep("LINE", names(my.varcorr), value=T)
    if(verbose) {
      cat("Genetic Variance = sum of", line.terms, "\n")
      cat("Environmental Variance = sum of", setdiff(names(my.varcorr), line.terms), "\n")
    }
    my.results["H2"] <- sum(my.varcorr[line.terms]) / sum(my.varcorr)
    
    my.rand.table <- rand(my.lmer)$rand.table
    pval.names <- paste0(sub(":",".",row.names(my.rand.table)), ".PVal")
    my.results[pval.names] <- my.rand.table[,"p.value"]
    
    # Also return the terms need to compute Coefficient of Genetic Variation
    # Which are the among line variance (just the numerator of H2 computation)
    # And the response mean (intercept of the model)
    # These two values need to be un-scaled before computing CGV, so just pass back here
    my.results["ALV"] <- sum(my.varcorr[line.terms])
    my.results["Mean"] <- fixef(my.lmer)["(Intercept)"]
  }
  names(my.results) <- gsub("[:]",".",names(my.results))
  
  # Compute line means
  # In all cases, EXCLUDE the Wolbachia effect
  
  # ---- MAJOR FIX - 2/20/18 ---- #
  # THIS IS WHERE THERE WAS AN ISSUE WITH INCLUDING CORRECTION FACTORS
  # THE ORIGINAL CODE WAS:
  #  remove.effects <- c("WOL","SEX:WOL",eqtl.inv, eqtl.pc)
  #  line.mean.terms <- setdiff(model.terms, remove.effects)
  #  line.mean.formula <- formula(paste0("~",paste(line.mean.terms, collapse="+")))
  #  if(verbose) {
  #    cat("Extracting line means using formula:\n")
  #    print(line.mean.formula)
  #  }
  #  my.fitted <- predict(my.lmer, re.form=line.mean.formula)
  # The intent was to drop the corrective factors out of the model formula
  # and then invoke predict to generate output files with the new reduced formula
  # HOWEVER, predict doesn't honor the re.form parameter, basically it always keeps all fixed effects from the original model
  # Wen suggested a simpler alternative, which is the following:
  # For within-sex baseline models (no treatment), just use ranef$LINE
  # For pooled-sex baseline models, ALSO use ranef$LINE
  #  This means we NEVER compute line means for each sex, we just take the sex-independent line means from the model
  # For within-sex treatment models (EtOH, 3WK), use ranef$TRMT:LINE
  #  This provides what we had before, which is a line mean for each treatment group
  #  We can still do correlations, differences, and any other tests that gets done off of those
  # For pooled-sex baseline models, do the same as for within-sex, so again we are getting the "sex-independent" portion of the model only
  # These should be very similar to what was initially intended, because the other terms would just cancel out or be fixed constants added to all line means
  if("TRMT" %in% my.model.terms) {
    my.line.means <- ranef(my.lmer, drop=T)$'TRMT:LINE'
    # The labels for these will be TRMT:LINE, but TRMT will be 0/1
    # They are NOT in same order my.model.input
    # We want terms to be in one of these formats:
    # LINE_TRMT (Pooled Sex)
    # LINE_SEX_TRMT (Within Sex)
    # This is mainly to keep the expected formatting for downstream functions
    # and to keep things matching the auto-generated GROUP column my.model.terms
    
    # Create a mapping of TRMT:LINE column values back to GROUP names in my.model.input
    term.group.table <- data.frame(
      MODEL=paste(my.model.input$TRMT, my.model.input$LINE, sep=":"),
      GROUP=my.model.input$GROUP
    )
    # For Pooled Sex models, drop the _F/M_ part of group names
    if("SEX" %in% my.model.terms) {
      term.group.table$GROUP <- sub("[_][FM]","",term.group.table$GROUP)
    }
    term.group.table <- unique(term.group.table)
    stopifnot(sum(duplicated(term.group.table$MODEL))==0)
    stopifnot(sum(duplicated(term.group.table$GROUP))==0)
    stopifnot(all(names(my.line.means) %in% term.group.table$MODEL))
    # Turn into a mapping from names(my.line.means) to groups
    term.group.map <- term.group.table$GROUP
    names(term.group.map) <- term.group.table$MODEL
    names(my.line.means) <- term.group.map[names(my.line.means)]
  } else {
    # In all other cases (including Pooled Sex), just extract the LINE terms (ignore Sex effects/Interactions)
    my.line.means <- ranef(my.lmer, drop=T)$LINE
    # IF This is a within-sex model then map LINE IDs to GROUP IDs in my.model.input
    # This keeps the column headers the same as the old code for downstream things
    if(!("SEX" %in% my.model.terms)) {
      line.group.table <- unique(my.model.input[,c("LINE","GROUP")])
      line.group.map <- line.group.table$GROUP
      names(line.group.map) <- line.group.table$LINE
      stopifnot(length(line.group.map)==length(my.line.means))
      stopifnot(all(names(line.group.map) %in% names(my.line.means)))
      names(my.line.means) <- line.group.map[names(my.line.means)]
    }
  }
  
  # Create a list structure that maps Groups to Samples, with groups ordered by line ID (and further grouped by SEX and/or TRMT if needed)
  if("TRMT" %in% colnames(my.model.input)) {
    if("SEX" %in% colnames(my.model.input)) {
      group.order <- order(my.model.input$TRMT, my.model.input$SEX, as.integer(my.model.input$LINE), decreasing=F)
    } else {
      group.order <- order(my.model.input$TRMT, as.integer(my.model.input$LINE), decreasing=F)
    }
  } else {
    if("SEX" %in% colnames(my.model.input)) {
      group.order <- order(my.model.input$SEX, as.integer(my.model.input$LINE), decreasing=F)
    } else {
      group.order <- order(as.integer(my.model.input$LINE), decreasing=F)
    }
  }
  
  my.groups <- unique(my.model.input$GROUP[group.order])
  # UPDATED 2/21/18 - Establishing group order is still useful here,
  # BUT we don't need to map group IDs to samples, we already have line means matched up to group IDs from the new code above
  # For pooled sex models, remove the _SEX portion of group names and collapse
  if("SEX" %in% my.model.terms) {
    my.groups <- sub("[_][MF]", "", my.groups)
    my.groups <- unique(my.groups)
  }
  stopifnot(length(my.groups)==length(my.line.means))
  stopifnot(all(my.groups %in% names(my.line.means)))
  my.line.means <- my.line.means[my.groups]
  # 2/21/18 - NO LONGER NEEDED
  # names(my.groups) <- my.groups
  # my.groups <- lapply(my.groups, function(x){row.names(my.model.input)[my.model.input$GROUP==x]})
  #
  # # Loop over groups, make sure all fitted values match, then collapse into table
  # my.line.means <- unlist(lapply(my.groups, function(x){
  #   if(length(x) == 1) {
  #     return(my.fitted[x])
  #   } else {
  #     # DON'T APPLY THIS CHECK DURING PERMUTATION TEST
  #     if(is.na(my.args$PERM)) {
  #       stopifnot(all(my.fitted[x] == my.fitted[x[1]]))
  #     }
  #     return(mean(my.fitted[x]))
  #   }
  # }))
  # names(my.line.means) <- names(my.groups)
  
  # For Treated models only, perform Levene's test for difference in variance between Treated and Untreated line means
  # --- DISABLED --- 
  # Now that line means come ONLY from LINE:TRMT,
  # this test no longer makes sense, and technically I don't think it should have been done here
  # The correct way to do this would be to run the within treatment models
  # (3wk and 5d models separately for same lines)
  # Then perform Levene's test on the line means from each of those models
  # if("TRMT" %in% my.model.terms) {
  #  ctrl.groups <- unique(my.model.input$GROUP[my.model.input$TRMT==0])
  #  trmt.groups <- unique(my.model.input$GROUP[my.model.input$TRMT==1])
  #  # For pooled sex model, drop sex component of GROUP names
  #  if("SEX" %in% my.model.terms) {
  #    ctrl.groups <- unique(sub("[_][FM]","",ctrl.groups))
  #    trmt.groups <- unique(sub("[_][FM]","",trmt.groups))
  #  }
  #  stopifnot(all(ctrl.groups %in% names(my.line.means)))
  #  stopifnot(all(trmt.groups %in% names(my.line.means)))
  #  stopifnot(length(intersect(ctrl.groups, trmt.groups))==0)
  #  # Perform Levene's test using function in car package
  #  my.results["Levene.PVal"] <- leveneTest(y=my.line.means[c(ctrl.groups,trmt.groups)], group=factor(c(rep("Ctrl", times=length(ctrl.groups)), rep("Trmt", times=length(trmt.groups)))))["group","Pr(>F)"]
  # }
  
  # RETURN LIST WITH MEMBERS model.results, line.means
  return(list(model.results=my.results, line.means=my.line.means))
}

# Function that takes a list of results from modelGeneticVar and compiles the model result or line mean table
# List member names should be the Gene/Feature IDs
# Set extract=line.means to get line means table, default is model.results table
extractModelTable <- function(model.outputs, extract="model.results") {
  my.table <- foreach(i=1:length(model.outputs), .combine='rbind') %do% {model.outputs[[i]][[extract]]}
  # If my.table has one row, it will come back as a vector, in which case it needs to be forced into being a data frame
  if(is.null(dim(my.table))) {
    my.table <- as.data.frame(t(my.table))
  }
  # Copy over row (feature) names, if present
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
    pval.cols <- grep("[.][Pp][Vv]al", colnames(my.table), value=T)
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
  for(j.name in pval.cols) {
    # Add corrected version at end of table
    # Hits an error if the column already exists!
    j.adj.col <- paste0(j.name, suffix)
    stopifnot(!(j.adj.col %in% colnames(my.table)))
    my.table[,j.adj.col] <- p.adjust(my.table[,j.name], method=method)
    
    # Get col number of original p-value column
    j <- which(colnames(my.table)==j.name)[1]
    # Adjusted column is always the last one
    # Reorder to put new column right after j
    if(j < (ncol(my.table)-1)) my.table <- my.table[,c(1:j,ncol(my.table),(j+1):(ncol(my.table)-1))]
  }
  
  # Return the expanded table
  return(my.table)
}



# --- Primary Model/Input Setup --- #

# This block sets up the terms used in the primary model and the subsets of samples to run the model on
# This block also determines the group names that are used in line mean output table

# All models have (1|LINE) random effect term
model.terms <- c("(1|LINE)")

# Decide whether to include WOL model term
if(my.args$WOLADJ) {
  model.terms <- c(model.terms, "WOL")
}

# For eQTL mode, include Inv and PC terms
if(my.args$EQTL) {
  model.terms <- c(model.terms, eqtl.inv, eqtl.pc)
}

# Decided whether to include TRMT and TRMT:LINE ixn terms
if(use.trmt) {
  model.terms <- c(model.terms, "TRMT", "(1|TRMT:LINE)")
}

# TO DO: Should there be a SEX:TRMT interaction term? Probably...

# If SEX column exists in sample.table, split the expr table by sex
sample.subsets <- list()
if(my.args$POOLED) {
  if("SEX" %in% colnames(sample.table)) {
    # Include SEX fixed effect and SEX:LINE random ixn effect
    model.terms <- c(model.terms, "SEX", "(1|SEX:LINE)")
    # Also include SEX:WOL term if doing Wolbachia correction
    if(my.args$WOLADJ) {
      model.terms <- c(model.terms, "SEX:WOL")
    }
    
    # Determine which lines have at least 2 replicates in EACH sex
    # this determines the one subset of samples that will be processed
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
    
    # If modelling trmt effects, include that in group names
    # Also include SEX:TRMT (fixed ixn) and SEX:TRMT:LINE 3-way random ixn terms here if needed
    if(use.trmt) {
      sample.table$GROUP <- paste(sample.table$LINE, sample.table$SEX, sample.table[,my.args$TRMT], sep="_")
      model.terms <- c(model.terms, c("SEX:TRMT", "(1|SEX:TRMT:LINE)"))
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

# Validate model and enforce term order
all.terms <- c("SEX","WOL","SEX:WOL",eqtl.inv,eqtl.pc,"TRMT","SEX:TRMT","(1|LINE)","(1|SEX:LINE)","(1|TRMT:LINE)","(1|SEX:TRMT:LINE)")
stopifnot(all(model.terms %in% all.terms))
model.terms <- intersect(all.terms, model.terms)
cat("\nFitting the following model for each expression feature:\n")
cat(getModelString(model.terms), "\n")
if(my.args$EQTL & my.args$POOLED) {
  cat("WARNING: Pooled eQTL model does not currently include SEX by INV/PC interaction terms!\n")
}
input.terms <- getInputTerms(model.terms)
stopifnot(all(input.terms %in% colnames(sample.table)))
cat("Required terms are present:", input.terms, "\n\n")


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
  # DEPRECATING THIS FILTER - SHOULD USE FLAGS AS FILTER INSTEAD
  # if(my.args$TYPE != "rle") {
  #  min.rep.expr.perc <- apply(ss.expr, 1, function(x){mean(x > log2(my.args$PSEUDO))})
  #  non.rep.features <- names(min.rep.expr.perc)[min.rep.expr.perc < 0.05]
  #  
  #  cat("Removing", length(non.rep.features), "features with expression in <",my.args$FDR*100,"% of this subset of samples.\n")
  #  keep.features <- setdiff(row.names(ss.expr), non.rep.features)
  #  ss.expr <- ss.expr[keep.features,]
  # }
  # cat("\n")
  
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
  # If linereg.table is present, do simple linear regression on that first, and take the residuals
  if(is.null(linereg.table)) {
    ss.linereg <- NULL
  } else {
    # Just transform linereg.table to have same dimensions as ss.expr here
    # Actual modeling of effects will be handled in core modeling function

    # Construct a vector that maps each sample to the appropriate column in linereg.table
    ss.samp.linereg.col <- paste0("line_", ss.samples$LINE)
    names(ss.samp.linereg.col) <- row.names(ss.samples)
    if(!all(ss.samp.linereg.col %in% colnames(linereg.table))) {
      stop("Missing line regression data for one or more sample lines.")
    }
    stopifnot(all(names(ss.samp.linereg.col)==row.names(ss.expr)))
    # Create variant rate table with same dimensions as ss.expr
    #  - just duplicate line values for replicate samples
    ss.linereg <- linereg.table[colnames(ss.expr),ss.samp.linereg.col]
    colnames(ss.linereg) <- names(ss.samp.linereg.col)
    ss.linereg <- t(ss.linereg)
    # Apply scaling here
    ss.linereg <- scale(ss.linereg)
    # Replace NA columns with all 0
    null.cols <- which(apply(ss.linereg, 2, function(x){all(is.na(x))}))
    if(length(null.cols) > 0) {
      ss.linereg[,null.cols] <- 0
    }
    stopifnot(sum(is.na(ss.linereg))==0)
    stopifnot(all(dim(ss.expr)==dim(ss.linereg)))
    stopifnot(all(row.names(ss.expr)==row.names(ss.linereg)))
    stopifnot(all(colnames(ss.expr)==colnames(ss.linereg)))
  }
  
  # 11/30/16 - Moved this down to be done AFTER the linereg step, so that it does NOT disrupt LINEREG correction
  # When PERM optin specified, 
  # Permute the sample IDs to get bg distribution of H2 values and p-values
  # TO DO: THIS IS NOT VALID FOR POOLED MODEL!!!!
  #   A better way to do this would be to define a bunch of possible permutations in a separate script and output to a table
  #   Then load that table here and apply a single permutation
  samp.shuf <- NULL
  if(!is.na(my.args$PERM)) {
    if(my.args$POOLED) {
      stop("PERM method is not yet set up for POOLED model!")
    }
    # Create a random ordering
    cat("Shuffling LINE IDs...\n")
    samp.shuf <- sample.int(n=nrow(ss.expr), replace=F)
    real.groups <- lapply(unique(ss.samples$LINE), function(x){which(ss.samples$LINE==x)})
    perm.samples <- ss.samples[rownames(ss.expr)[samp.shuf],]
    perm.groups <- lapply(unique(perm.samples$LINE), function(x){which(perm.samples$LINE==x)})
    # Compute how many biological replicates are still paired together after shuffling
    # Should be 0...
    unperm.groups <- sum(unlist(lapply(perm.groups, function(x){
      max(unlist(lapply(real.groups, function(y){length(intersect(x,y))}))) - 1
    })))
    if(unperm.groups > 0) {
      cat("WARNING:", unperm.groups, "true replicate pairs are still paired after shuffling.\n\n")
    }
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
    i.model.input <- getModelInput(i=i, my.expr=ss.expr, my.samples=ss.samples, line.reg=ss.linereg, model.inputs=input.terms)
    # Scramble JUST the LINE/GROUP terms
    if(!is.null(samp.shuf)) {
      i.model.input$LINE <- i.model.input$LINE[samp.shuf]
      i.model.input$GROUP <- i.model.input$GROUP[samp.shuf]
    }
    modelGeneticVar(i.model.input, model.terms)
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
  effect.cols <- c(grep("[.]Coeff$", colnames(ss.model.results), value=T), grep("[.]Log2FC$", colnames(ss.model.results), value=T))
  for(fc.col in effect.cols) {
    ss.model.results[,fc.col] <- ss.model.results[,fc.col] * ss.sd
  }
  
  # For ALV and Mean, the unscaling is a bit more complicated...
  # ALV = Among Line Variance, which was the numerator in H^2 calculation
  # Since this is a variance measure, it should be multiplies by ss.sd ** 2
  if("ALV" %in% colnames(ss.model.results)) {
    ss.model.results[,"ALV"] <- ss.model.results[,"ALV"] * (ss.sd ** 2)
  }
  # Mean gets unscaled just like line means
  if("Mean" %in% colnames(ss.model.results)) {
    ss.model.results[,"Mean"] <- (ss.model.results[,"Mean"] * ss.sd) + ss.center
  }
  
  # Now compute the Coefficient of Genetic Variation (CGV)
  # CGV = (100*sqrt(ALV))/Mean
  # See: https://ecoqui.wordpress.com/2015/03/06/how-to-calculate-genetic-variance-components-coefficient-of-genetic-variation-and-genetic-correlations-in-r/
  # And: Felix, et al. 2012
  if(all(c("ALV","Mean") %in% colnames(ss.model.results))) {
    ss.model.results[,"CGV"] <- 100*sqrt(ss.model.results[,"ALV"])/ss.model.results[,"Mean"]
    # Drop ALV and Mean here
    ss.model.results <- ss.model.results[,setdiff(colnames(ss.model.results),c("ALV","Mean"))]
  }
  
  # Perform multiple testing correction on all P-value columns
  ss.model.results <- correctPvalCols(ss.model.results, method="BH")
  
  # Report LINEREG correction parameters if computed
  if("LINEREG.Coeff" %in% colnames(ss.model.results)) {
    cat(sum(ss.model.results$LINEREG.PVal.FDR <= my.args$FDR),"features have significant Line Regression effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
    cat(" =", round(mean(ss.model.results$LINEREG.PVal.FDR <= my.args$FDR)*100), "% of all features tested.\n\n")
    cat("Distribution of % Variance explained by Line Regression data:\n")
    cat("Range:", range(ss.model.results[,"LINEREG.Perc.Var"]), "\n")
    cat("Box Stats:", boxplot.stats(ss.model.results[,"LINEREG.Perc.Var"])$stats,"\n")
    # Looking at raw log2FC is really meaningless, scrubbing this from output...
    # cat("Distribution of Line Regression Effect Sizes (Log2FC):\n")
    # cat("Range:", range(ss.model.results[,"LINEREG.Log2FC"]), "\n")
    # cat("Box Stats:", boxplot.stats(ss.model.results[,"LINEREG.Log2FC"])$stats,"\n")
    cat("\n")
  }
  
  # Report significant Sex effects:
  if(my.args$POOLED) {
    signif.rows <- ss.model.results[,"SEX.PVal.FDR"] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    signif.rows[failed.rows] <- F
    if(sum(failed.rows) > 0) {
      cat(sum(failed.rows),"features failed to converge when testing Sex fixed effect ANOVA\n")
      cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
    }
    cat(sum(signif.rows),"features have significant Sex effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
    cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n\n")
    # cat("Distribution of Significant Sex Effects (Log2):\n")
    # cat("Range:", range(ss.model.results[signif.rows,"SEX.Log2FC"]), "\n")
    # cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"SEX.Log2FC"])$stats,"\n")
    cat("Distribution of % Variance Explained by Significant Sex Effects:\n")
    cat("Range:", range(ss.model.results[signif.rows,"SEX.Perc.Var"]), "\n")
    cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"SEX.Perc.Var"])$stats,"\n")
    cat("\n")
  }
  
  # Report significant Wolbachia effects:
  if(my.args$WOLADJ) {
    # Consider WOL effect significance (and SEX:WOL effect for pooled model)
    signif.rows <- ss.model.results[,"WOL.PVal.FDR"] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    signif.rows[failed.rows] <- F
    if(sum(failed.rows) > 0) {
      cat(sum(failed.rows),"features failed to converge when testing Wolbachia fixed effect ANOVA\n")
      cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
    }
    cat(sum(signif.rows),"features have significant Wolbachia effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
    cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n\n")
    if(any(signif.rows)) {
      # cat("Distribution of Significant Wolbachia Effects (Log2):\n")
      # cat("Range:", range(ss.model.results[signif.rows,"WOL.Log2FC"]), "\n")
      # cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"WOL.Log2FC"])$stats,"\n")
      cat("Distribution of % Variance Explained by Significant Wolbachia Effects:\n")
      cat("Range:", range(ss.model.results[signif.rows,"WOL.Perc.Var"]), "\n")
      cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"WOL.Perc.Var"])$stats,"\n")
    }
    cat("\n")
      
    if(my.args$POOLED) {
      signif.rows <- ss.model.results[,"SEX.WOL.PVal.FDR"] <= my.args$FDR
      failed.rows <- is.na(signif.rows)
      signif.rows[failed.rows] <- F
      if(sum(failed.rows) > 0) {
        cat(sum(failed.rows),"features failed to converge when testing Sex:Wolbachia fixed interaction effect ANOVA\n")
        cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
      }
      cat(sum(signif.rows),"features have significant Sex:Wolbachia effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
      cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n\n")
      if(any(signif.rows)) {
        # cat("Distribution of Sex:Wolbachia interaction Effects (Log2):\n")
        # cat("Range:", range(ss.model.results[signif.rows,"SEX.WOL.Log2FC"]), "\n")
        # cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"SEX.WOL.Log2FC"])$stats,"\n")
        cat("Distribution of Significant Sex:Wolbachia % Variance Explained:\n")
        cat("Range:", range(ss.model.results[signif.rows,"SEX.WOL.Perc.Var"]), "\n")
        cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"SEX.WOL.Perc.Var"])$stats,"\n")
      }
      cat("\n")
    }
  }
  
  if(my.args$EQTL) {
    # Check for significant effects for each Inversion and PC term
    for(j in c(eqtl.inv, eqtl.pc)) {
      j.fdr.col <- paste0(j, ".PVal.FDR")
      signif.rows <- ss.model.results[,j.fdr.col] <= my.args$FDR
      failed.rows <- is.na(signif.rows)
      signif.rows[failed.rows] <- F
      if(sum(failed.rows) > 0) {
        cat(sum(failed.rows),"features failed to converge when testing", j, "fixed effect ANOVA\n")
        cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
      }
      cat(sum(signif.rows),"features have significant", j, "effects for", ss, "samples at", my.args$FDR*100,"% FDR.\n")
      cat(" =", round(mean(signif.rows)*100),"% of all features tested.\n\n")
      if(any(signif.rows)) {
        j.pv.col <- paste0(j, ".Perc.Var")
        cat("Distribution of Significant",j,"% Variance Explained:\n")
        cat("Range:", range(ss.model.results[signif.rows,j.pv.col]), "\n")
        cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,j.pv.col])$stats, "\n")
      }
      cat("\n")
    }
  }
  
  if(use.trmt) {
    fdr.col <- "TRMT.PVal.FDR"
    fc.col <- "TRMT.Log2FC"
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
  
  if(use.trmt & my.args$POOLED) {
    fdr.col <- "SEX.TRMT.PVal.FDR"
    signif.rows <- ss.model.results[,fdr.col] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    signif.rows[failed.rows] <- F
    if(sum(failed.rows) > 0) {
      cat(sum(failed.rows),"features failed to converge when testing Sex by Treatment fixed interaction ANOVA\n")
      cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
    }
    cat(sum(signif.rows),"features have significant Sex:Treatment interaction effects for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
    cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n\n")
    cat("\n")
  }
  
  # Report significant genetic variance 
  fdr.cols <- c("LINE.PVal.FDR", "SEX.LINE.PVal.FDR", "TRMT.LINE.PVal.FDR", "SEX.TRMT.LINE.PVal.FDR")
  fdr.cols <- intersect(fdr.cols, colnames(ss.model.results))
  stopifnot(length(fdr.cols) > 0)
  if(length(fdr.cols) == 1) {
    signif.rows <- ss.model.results[,fdr.cols] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    signif.rows[failed.rows] <- F
  } else {
    signif.rows <- apply(ss.model.results[,fdr.cols], 1, function(x){any(x <= my.args$FDR, na.rm=F)})
    failed.rows <- apply(ss.model.results[,fdr.cols], 1, function(x){any(is.na(x))})
  }
  if(sum(failed.rows) > 0) {
    cat(sum(failed.rows),"feature models failed to converge for", ss, "samples\n")
    cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
  }
  cat(sum(signif.rows),"features have significant genetic variance in mixed effects model for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
  cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n")
  cat("Distribution of H^2 values for significantly heritable features:\n")
  cat("Range:", range(ss.model.results[signif.rows,"H2"]), "\n")
  cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"H2"])$stats,"\n")
  cat("\n")
  
  # If there's a TRMT effect, establish list of corresponding groups for treatment and controls
  if(use.trmt) {
    # NEW 2/21/18 - Matching up treated and untreated line means changed, depending on within sex vs pooled sex model
    if(ss == "Pooled") {
      # For Pooled Sex models, the line means are now already collapsed by sex
      # So drop the _F and _M component of group names and SEX column
      group.map.table <- ss.samples[,c("LINE","TRMT","GROUP")]
      group.map.table$GROUP <- sub("[_][MF]", "", group.map.table$GROUP)
      control.group.table <- unique(group.map.table[group.map.table$TRMT==0,])
      control.groups <- control.group.table$GROUP
      names(control.groups) <- control.group.table$LINE
      treated.group.table <- unique(group.map.table[group.map.table$TRMT==1,])
      treated.groups <- treated.group.table$GROUP
      names(treated.groups) <- treated.group.table$LINE
    } else {
      # The old method still applies for within sex models
      control.groups <- unique(ss.samples[ss.samples$TRMT==0,"GROUP"])
      control.group.lines <- unlist(lapply(control.groups, function(x){unique(ss.samples[ss.samples$GROUP==x,"LINE"])}))
      control.group.sex <- unlist(lapply(control.groups, function(x){unique(ss.samples[ss.samples$GROUP==x,"SEX"])}))
      control.group.line.sex <- paste(control.group.lines, control.group.sex, sep="_")
      names(control.groups) <- control.group.line.sex
      control.groups <- control.groups[sort(names(control.groups))]
    
      treated.groups <- unique(ss.samples[ss.samples$TRMT==1,"GROUP"])
      treated.group.lines <- unlist(lapply(treated.groups, function(x){unique(ss.samples[ss.samples$GROUP==x,"LINE"])}))
      treated.group.sex <- unlist(lapply(treated.groups, function(x){unique(ss.samples[ss.samples$GROUP==x,"SEX"])}))
      treated.group.line.sex <- paste(treated.group.lines, treated.group.sex, sep="_")
      names(treated.groups) <- treated.group.line.sex
      treated.groups <- treated.groups[sort(names(treated.groups))]
    }
    
    # TEMP DEBUG
    # cat("control.groups before sorting:\n")
    # print(head(control.groups))
    # cat("treated.groups before sorting:\n")
    # print(head(treated.groups))
    
    # Make sure both are sorted the same and then match up
    control.groups <- control.groups[order(names(control.groups))]
    treated.groups <- treated.groups[order(names(treated.groups))]
    
    # TEMP DEBUG
    # cat("control.groups after soring:\n")
    # print(head(control.groups))
    # cat("treated.groups after sorting:\n")
    # print(head(treated.groups))
    
    stopifnot(sum(duplicated(control.groups))==0)
    stopifnot(sum(duplicated(names(control.groups)))==0)
    stopifnot(sum(duplicated(treated.groups))==0)
    stopifnot(sum(duplicated(names(treated.groups)))==0)
    stopifnot(length(control.groups)==length(treated.groups))
    stopifnot(all(names(control.groups)==names(treated.groups)))
    stopifnot(all(control.groups %in% colnames(ss.line.means)))
    stopifnot(all(treated.groups %in% colnames(ss.line.means)))
    stopifnot(length(intersect(control.groups,treated.groups))==0)
  
    # --- DISABLED --- 
    # Now that line means = TRMT:LINE terms only, correlations do not give the intended result
    # Technically, the better way to do this is run the corresponding within treatment models
    # Then compute correlations of the line means that come out of those
    #   
  #   # Compute treated vs control line mean correlations for each gene
  #   # Can use parallelization here
  #   ss.cor <- foreach(i=1:nrow(ss.line.means), .combine='rbind') %dopar% {
  #     # Just computing PCC for now, but can easily add SCC here
  #     # Also reporting the correlation p-value
  #     i.cor.test <- cor.test(t(ss.line.means[i,treated.groups]), t(ss.line.means[i,control.groups]), method="pearson")
  #     return(data.frame(PCC=i.cor.test$estimate, PCC.PVal=i.cor.test$p.value))
  #   }
  #   row.names(ss.cor) <- row.names(ss.line.means)
  #   
  #   # Benjamini Hochberg correction of Correlation P-values
  #   ss.cor <- correctPvalCols(ss.cor, method="BH")
  #   ss.model.results <- cbind(ss.model.results, ss.cor)
  #   
  #   # Report some overall stats here
  #   fdr.col <- "PCC.PVal.FDR"
  #   signif.rows <- ss.model.results[,fdr.col] <= my.args$FDR
  #   failed.rows <- is.na(signif.rows)
  #   signif.rows[failed.rows] <- F
  #   if(sum(failed.rows) > 0) {
  #     cat(sum(failed.rows),"features have invalid line means for computing PCC.\n")
  #     cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
  #   }
  #   cat(sum(signif.rows),"features have significant PCC between Treated and Control groups for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
  #   cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n")
  #   cat("Distribution of PCC values for significantly correlated features:\n")
  #   cat("Range:", range(ss.model.results[signif.rows,"PCC"]), "\n")
  #   cat("Box Stats:", boxplot.stats(ss.model.results[signif.rows,"PCC"])$stats,"\n\n")
  }
  
  # If Levene's test was performed, report the number of features passing FDR
  if("Levene.PVal.FDR" %in% colnames(ss.model.results)) {
    signif.rows <- ss.model.results[,"Levene.PVal.FDR"] <= my.args$FDR
    failed.rows <- is.na(signif.rows)
    if(sum(failed.rows) > 0) {
      cat(sum(failed.rows),"features have invalid line means for performing Levene's test.\n")
      cat(" =", round(mean(failed.rows)*100), "% of all features tested.\n")
    }
    cat(sum(signif.rows),"features have significant difference in variance between Treated and Control groups for", ss, "samples at",my.args$FDR*100,"% FDR.\n")
    cat(" =", round(mean(signif.rows)*100), "% of all features tested.\n\n")
  }
  
  # If CGV was computed, just give the overall distribution
  # TO DO: Maybe should only give this for those with significantly heritable features?
  if("CGV" %in% colnames(ss.model.results)) {
    cat("Distribution of CGV values for all features:\n")
    cat("Range:", range(ss.model.results[,"CGV"]), "\n")
    cat("Box Stats:", boxplot.stats(ss.model.results[,"CGV"])$stats,"\n\n")
  }
    
  # Store this table in ss.results list
  ss.results[[ss]] <- ss.model.results
  
  # Output line means ONLY for non-permutation analysis
  if(my.args$MEANS) {
    # UNSCALE HERE
    for(j in 1:ncol(ss.line.means)) {
      ss.line.means[,j] <- (ss.line.means[,j] * ss.sd) + ss.center
    }
    
    # Create output table, first column is gene ID, each subsequent column is group mean expr fitted from regression
    line.mean.table <- cbind(GENE=row.names(ss.line.means), FLAG=gene.flags[row.names(ss.line.means)], ss.line.means)
    
    # Report ONLY features with significant LINE term here
    # This should not impact the line diffs below,
    # which will output all features with TRMT:LINE ixn, regardless of LINE p-value
    # Also, because sex effects/ixns are now always ignored for line means from pooled models
    # there is no reason to consider those interactions here
    # test.cols <- paste0(c("LINE", "SEX.LINE", "TRMT.LINE", "SEX.TRMT.LINE"), ".PVal.FDR")
    test.cols <- paste0(c("LINE"), ".PVal.FDR")
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
    line.mean.file <- paste0(my.args$OUTSTUB, "_", ss, "_line_means.txt")
    cat("Writing line means for", ss, "samples to:", line.mean.file, "\n")
    if(file.exists(line.mean.file)) {
      cat("WARNING: Overwriting existing file\n")
    }
    write.table(line.mean.table, line.mean.file, sep="\t", row.names=F, quote=F)
    
    # If use.trmt=TRUE, also compute another table with line differences (Treated-Control)
    if(use.trmt) {
      line.diff.table <- ss.line.means[,treated.groups]-ss.line.means[,control.groups]
      colnames(line.diff.table) <- names(treated.groups)
      line.diff.table <- cbind(GENE=row.names(ss.line.means), FLAG=gene.flags[row.names(ss.line.means)], line.diff.table)
      
      # Further filter this table to just those with significant TRMT.LINE interactions
      # Remove genes that were not significantly heritable here
      # Originally this would output line diffs for genes with TRMT:LINE or SEX:TRMT:LINE interactions
      # But because no SEX effects/interactions are included in line mean computation anymore
      # there is no reason to include features based on SEX:TRMT:LINE interaction here
      test.cols <- paste0(c("TRMT.LINE"), ".PVal.FDR")
      test.cols <- intersect(test.cols, colnames(ss.model.results))
      stopifnot(length(test.cols) > 0)
      
      if(length(test.cols) == 1) {
        line.diff.rows <- row.names(ss.model.results)[ss.model.results[,test.cols] <= my.args$FDR]
      } else {
        line.diff.rows <- row.names(ss.model.results)[apply(ss.model.results[,test.cols], 1, function(q){any(q <= my.args$FDR)})]
      }
      cat("Outputting line diffs for", length(line.diff.rows), "features with FDR <", my.args$FDR, "in any of these columns:\n")
      cat(test.cols, "\n")
      
      line.diff.table <- line.diff.table[line.diff.rows,]
      
      # Output the table to disk, overwrite if it exists
      line.diff.file <- paste0(my.args$OUTSTUB, "_", ss, "_line_diffs.txt")
      cat("Writing line diffs for", ss, "samples to:", line.diff.file, "\n")
      if(file.exists(line.diff.file)) {
        cat("WARNING: Overwriting existing file\n")
      }
      write.table(line.diff.table, line.diff.file, sep="\t", row.names=F, quote=F)
      
    }
    
    # UPDATED - 2/21/18
    # This part is no longer needed
    # The new method for extracting line means from Pooled Sex models is already independent of sex.
    # So now the Pooled_line_means.txt file should be used for downstream analysis that used to use Pooled_SexAvg_line_means.txt
    # # If POOLED=TRUE, also compute another table with average line means across sexes (for features with no sex-line interactions, this is the type of line mean to use)
    # if(my.args$POOLED) {
    #   if(!use.trmt) {
    #     # When not treated samples, just group by line
    #     # Create list structure of all groups (column names in line.mean.table) mapping to each LINE
    #     all.lines <- unique(ss.samples$LINE)
    #     line.groups <- lapply(all.lines, function(line){unique(ss.samples$GROUP[ss.samples$LINE==line])})
    #     names(line.groups) <- all.lines
    #     stopifnot(all(unlist(line.groups) %in% colnames(line.mean.table)))
    #     stopifnot(all(unlist(lapply(line.groups, length)) > 1))
    #     # NOTE: This keeps the same "significant" genes/feature set as for the individual line_sex mean table
    #     sex.avg.table <- line.mean.table[,c("GENE","FLAG")]
    #     for(line in names(line.groups)) {
    #       sex.avg.table[,line] <- apply(line.mean.table[,line.groups[[line]]], 1, mean)
    #     }
    #     sex.avg.file <- paste0(my.args$OUTSTUB, "_SexAvg_line_means.txt")
    #     cat("Writing Sex Avg line means to:", sex.avg.file, "\n")
    #     if(file.exists(sex.avg.file)) {
    #       cat("WARNING: Overwriting existing file\n")
    #     }
    #     write.table(sex.avg.table, sex.avg.file, sep="\t", row.names=F, quote=F)
    #   } else {
    #     # If treated, need to compute averages within line_trmt groups, then compute differences
    #     all.lines <- sort(unique(ss.samples$LINE))
    #     line.control.groups <- lapply(all.lines, function(line){unique(ss.samples$GROUP[(ss.samples$LINE==line)&(ss.samples$TRMT==0)])})
    #     names(line.control.groups) <- paste(all.lines, trmt.groups[1], sep="_")
    #     line.treated.groups <- lapply(all.lines, function(line){unique(ss.samples$GROUP[(ss.samples$LINE==line)&(ss.samples$TRMT==1)])})
    #     names(line.treated.groups) <- paste(all.lines, trmt.groups[2], sep="_")
    #     
    #     stopifnot(all(unlist(line.control.groups) %in% colnames(line.mean.table)))
    #     stopifnot(all(unlist(lapply(line.control.groups, length)) > 1))
    #     stopifnot(all(unlist(line.treated.groups) %in% colnames(line.mean.table)))
    #     stopifnot(all(unlist(lapply(line.treated.groups, length)) > 1))
    #     
    #     # NOTE: This keeps the same "significant" genes/feature set as for the individual line_sex mean table
    #     sex.avg.table <- line.mean.table[,c("GENE","FLAG")]
    #     for(line in names(line.control.groups)) {
    #       sex.avg.table[,line] <- apply(line.mean.table[,line.control.groups[[line]]], 1, mean)
    #     }
    #     for(line in names(line.treated.groups)) {
    #       sex.avg.table[,line] <- apply(line.mean.table[,line.treated.groups[[line]]], 1, mean)
    #     }
    #     sex.avg.file <- paste0(my.args$OUTSTUB, "_SexAvg_line_means.txt")
    #     cat("Writing Sex Avg line means to:", sex.avg.file, "\n")
    #     if(file.exists(sex.avg.file)) {
    #       cat("WARNING: Overwriting existing file\n")
    #     }
    #     write.table(sex.avg.table, sex.avg.file, sep="\t", row.names=F, quote=F)
    #     
    #     # Now compute the differences based off the averages
    #     sex.avg.diff.table <- sex.avg.table[,names(line.treated.groups)] - sex.avg.table[,names(line.control.groups)]
    #     colnames(sex.avg.diff.table) <- all.lines
    #     # Subset to genes in line.diff.table and add back GENE, FLAG cols
    #     sex.avg.diff.table <- cbind(line.diff.table[,c("GENE","FLAG")], sex.avg.diff.table[row.names(line.diff.table),])
    #     
    #     sex.avg.diff.file <- paste0(my.args$OUTSTUB, "_SexAvg_line_diffs.txt")
    #     cat("Writing Sex Avg line diffs to:", sex.avg.diff.file, "\n")
    #     if(file.exists(sex.avg.diff.file)) {
    #       cat("WARNING: Overwriting existing file\n")
    #     }
    #     write.table(sex.avg.diff.table, sex.avg.diff.file, sep="\t", row.names=F, quote=F)
    #   }
    # }
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
fdr.cols <- c(LINE="LINE.PVal.FDR")
if(my.args$POOLED) {
  fdr.cols <- c(fdr.cols, 'SEX:LINE'="SEX.LINE.PVal.FDR")
}
if(use.trmt) {
  fdr.cols <- c(fdr.cols, 'TRMT:LINE'="TRMT.LINE.PVal.FDR")
}
if(my.args$POOLED & use.trmt) {
  fdr.cols <- c(fdr.cols, 'SEX:TRMT:LINE'="SEX.TRMT.LINE.PVal.FDR")
}
fdr.col.names <- names(fdr.cols)

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
write.table(cbind(GENE=row.names(all.ss.results), FLAG=gene.flags[row.names(all.ss.results)], all.ss.results), my.args$H2FILE, sep="\t", row.names=F, quote=F)

cat("\n\nScript completed successfully!\n\n")
print(proc.time())
quit("no")
