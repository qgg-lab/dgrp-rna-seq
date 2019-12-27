#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 3/9/17
#
# trait_regression.R
#
# Goal: Load in expr feature line means and assess their predictive value for all available line features
# 
# Usage:
# trait_regression.R [OPTIONS] EXPR=subdir/..._line_means.txt
#  EXPR=  The file containing expression values to assess heritability on
#  SEX=   F/M (Default=infer from EXPR file), determines which sex to compare to for line means
#  TRAITDIR=   Path to trait file (Defaults to /home/ljeveret/Resources/DGRP/Freeze2/)
#  TRAITFILE=  Name of trait data file for appropiate sex, 
#   defaults to ~/Resources/DGRP/Freeze2/line_means_Wol_(female|male).txt
#  OUTDIR=  Path to store output data
#   defaults to the path to EXPR input, drop any trailing genVar/ or pca/, add traitreg/
#  OUTPUT=  Prefix for the output files, defaults to EXPR basename minus _line_means.txt (PCs left in to differentiate)
#
# Automatically puts all output in a traitCorr subdirectory in the same parent directory as data file, with similar basenames
#
# TO DO: Parallelize with doMC
# TO DO: Scale/center expression profiles for numerical accuracy?
# TO DO: Create a separate directory for each trait? - This might become important as more things are reported
# TO DO: Would also be worth doing multiple regression (all features? signif only? forward?) and report the total % of variance explained, 
#        plot predicted vs actual values too?
#
# TO DO: Option to only run against a subset of traits (or can just create alternate versions of the trait file with subset of cols)
# TO DO: Remove LOG2 option? Or make it an option to Log2 transform the traits only?
# TO DO: Create a parameter to make the Wolbachia adjusted version optional?
# TO DO: Should perform permutations to compute non-parametric FDR?
#

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

options(stringsAsFactors=F)

# -- Fixed Parameters -- #
# These could be made adjustable in the command-line args
pseudo.count <- 0.001
fdr.cutoff <- 0.05

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript trait_correlations.R [OPTIONS] EXPR=subdir/..._line_means.txt"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("EXPR=microbiome_species_filtered/pca/combined_microbe_filtered_rpm_WolAdj_SexAvg_PCs.txt", "TRAITDIR=~/Resources/DGRP/Freeze2/", "SEX=F")
if(length(my.args) == 0) {
  stop("Requires at least one argument. ", usageStr)
} else {
  cat("Running with params:",my.args,"\n")
}
# Parse param names and values
argSplit=strsplit(my.args, split="=", fixed=T)
my.args <- lapply(argSplit, "[", 2)
names(my.args) <- unlist(lapply(argSplit, "[", 1))

# Make sure EXPR param exists
if(!("EXPR" %in% names(my.args))) {
  stop("Missing EXPR parameter. ", usageStr)
}

# Make sure EXPR file exists
if(!file.exists(my.args$EXPR)) {
  stop("Specified EXPR file does not exist: ", my.args$EXPR)
}

# Fill in default value for missing params
setDefault <- function(param, default) {
  if(!(param %in% names(my.args))) {
    cat("Setting ",param,"=",default," by default\n", sep="")
    my.args[[param]] <<- default
  }
}

setDefault("SEX","auto")
if(my.args$SEX=="auto") {
  # Try to auto-determine SEX parameter from EXPR filename
  if(grepl("[_]F[_]line[_]means.txt$", my.args$EXPR) | grepl("[_]F[_]PCs[.]txt$", my.args$EXPR)) {
    my.args$SEX <- "F"
  } else if(grepl("[_]M[_]line[_]means.txt$", my.args$EXPR) | grepl("[_]M[_]PCs[.]txt$", my.args$EXPR)) {
    my.args$SEX <- "M"
  } else {
    cat("WARNING: Could not auto-determine SEX from EXPR file, using SEX=F by default.\n")
    my.args$SEX <- "F"
  }
}
# Make sure SEX=F or M
if(!(my.args$SEX %in% c("F","M"))) {
  cat("WARNING: Did not recognize SEX =", my.args$SEX, "using SEX=F by default.\n")
  my.args$SEX <- "F"
}

setDefault("TRAITDIR","/home/ljeveret/Resources/DGRP/Freeze2/")
if(!grepl("/$",my.args$TRAITDIR)) {
  my.args$TRAITDIR <- paste0(my.args$TRAITDIR, "/")
}

default.traitfile <- paste0("line_means_Wol_", if(my.args$SEX=="F"){"fe"}, "male.txt")
setDefault("TRAITFILE",default.traitfile)

# Append TRAITDIR to TRAITFILE, make sure it exists
my.args$TRAITFILE <- paste0(my.args$TRAITDIR, my.args$TRAITFILE)
if(!file.exists(my.args$TRAITFILE)) {
  stop("Specified TRAITFILE does not exist: ", my.args$TRAITFILE)
}

# DEPRECATED
# setDefault("LOG2",F)
# # If "0" or "1" map these to "FALSE" and "TRUE"
# if(is.character(my.args$LOG2) & (my.args$LOG2 == "0")){my.args$LOG2 <- F}
# if(is.character(my.args$LOG2) & (my.args$LOG2 == "1")){my.args$LOG2 <- T}
# # Convert any other type to logical
# if(!is.logical(my.args$LOG2)){my.args$LOG2 <- as.logical(my.args$LOG2)}
# # If conversion produced NA, just set to F
# if(is.na(my.args$LOG2)) {
#   cat("WARNING: Did not recognize value of LOG2 parameter, setting to F by default.\n")
#   my.args$LOG2 <- F
# }
# # UPDATE: 6/20/16
# # LOG2 transformation is no longer applied here, it is applied directly to the line_means file
# # Generate a big warning for now
# # TO DO: this option should be completely removed, along with pseudo-count option!
# if(my.args$LOG2) {
#   cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
#   cat("!           WARNING           !\n")
#   cat("! LOG2 OPTION IS DEPRECATED   !\n")
#   cat("! Are you sure the line means !\n")
#   cat("! are not already log scale?  !\n")
#   cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
# }
  
my.args$PSEUDO <- pseudo.count
my.args$FDR <- fdr.cutoff

# Extract path to EXPR input file, use it as the output dir
# (If EXPR data is under a genVar subdirectory, drop that from OUTDIR path)
default.outdir <- dirname(my.args$EXPR)
default.outdir <- sub("genVar[/]?$","", default.outdir)
default.outdir <- sub("pca[/]?$","",default.outdir)
if(default.outdir=="") {
  default.outdir <- "."
}
if(!grepl("/$",default.outdir)) {
  default.outdir <- paste0(default.outdir,"/")
}
default.outdir <- paste0(default.outdir, "traitreg/")
setDefault("OUTDIR",default.outdir)
if(!grepl("/$",my.args$OUTDIR)) {
  my.args$OUTDIR <- my.args$OUTDIR
}

# Base name of COUNTS file comes from shaving off the trailing _[FM]_line_means.txt and adding _SEX
default.output <- basename(my.args$EXPR)
default.output <- sub("[.]txt$", "", default.output)
default.output <- sub("[_]line[_]means$", "", default.output)
# Check if Sex is already denoted in file name, otherwise put on the end
if(!grepl(paste0("[_]",my.args$SEX,"[_]"),default.output)) {
  default.output <- paste0(default.output,"_",my.args$SEX)
}
setDefault("OUTPUT",default.output)


# --- LOAD INPUT --- #

# Load the line mean table
expr.table <- read.table(my.args$EXPR, header=T, sep="\t", row.names=1)
# Remove 'X' that R automatically puts in front of column names
colnames(expr.table) <- sub("^X", "", colnames(expr.table))
features <- row.names(expr.table)
cat("Loaded counts for",length(features),"features across",ncol(expr.table),"lines from",my.args$EXPR,"\n")

trait.table <- read.table(my.args$TRAITFILE, header=T, sep="\t", row.names=1)
# Remove 'X' that R automatically puts in front of column names
colnames(trait.table) <- sub("^X", "", colnames(trait.table))
# Transpose so that columns = lines
trait.table <- as.data.frame(t(trait.table))
cat("Loaded phenotype values for",nrow(trait.table),"traits across",ncol(trait.table),"lines from",my.args$TRAITS,"\n")

# Change column header format to match trait.table
colnames(expr.table) <- paste0("DGRP_", sub(paste0("_",my.args$SEX), "", colnames(expr.table)))
# Match up columns, drop lines that aren't shared
common.lines <- intersect(colnames(expr.table), colnames(trait.table))
cat(length(common.lines), "lines match up between expression and trait tables.\n")
if(length(setdiff(colnames(expr.table), colnames(trait.table))) > 0) {
  cat("Dropping",length(setdiff(colnames(expr.table), colnames(trait.table))),"lines from expression table due to lack of trait data.\n")
}
if(length(setdiff(colnames(trait.table), colnames(expr.table))) > 0) {
  cat("Dropping",length(setdiff(colnames(trait.table), colnames(expr.table))),"lines from trait table due to lack of expression data.\n")
}
# Make sure column order matches up even if we don't need to drop columns
expr.table <- expr.table[,common.lines]
trait.table <- trait.table[,common.lines]

# DEPRECATED
# if(my.args$LOG2) {
#  cat("Converting expression data to Log2 scale with pseudo-count =",my.args$PSEUDO,"\n")
#  expr.table <- log2(expr.table+my.args$PSEUDO)
# }
cat("\n")



# --- COMPUTE REGRESSIONS --- #

# This table will store trait correlations and p-values for each gene
# [DEPRECATED] Output data for each trait separately:
# expr.trait.reg <- data.frame(row.names=row.names(expr.table))

# loop over each trait and fit an ANOVA to get % variance explained and p-value
cat("Testing % Var Explained and ANOVA P-values for",nrow(expr.table),"expression features against", nrow(trait.table), "traits.\n")

# TO DO: Could use foreach construct here, but better to parallelize on the individual genes?
for(trait in row.names(trait.table)) {
  # Extract trait values to a vector, drop the ones that are NA
  trait.vals <- as.numeric(trait.table[trait,])
  names(trait.vals) <- colnames(trait.table)
  trait.vals <- trait.vals[!is.na(trait.vals)]
  cat("Computing regression models for",trait,"trait data from",length(trait.vals),"lines.\n")
  expr.subset <- expr.table[,names(trait.vals)]
  
  # Loop over the regression models, return 
  # TO DO: Use foreach and %dopar% to speed this up?
  expr.reg <- lapply(row.names(expr.subset), function(gene){
    trait.gene.lm.input <- as.data.frame(cbind(TRAIT=trait.vals, EXPR=as.numeric(expr.subset[gene,])))
    # Fit Linear Model over just the top PCs
    trait.gene.lm <- lm(TRAIT ~ EXPR, data=trait.gene.lm.input)
    trait.gene.coeff <- trait.gene.lm$coefficients["EXPR"]
    
    # Perform ANOVA to get P-value for each PC against this gene expression profile
    trait.gene.anova.pval <- as.data.frame(anova(trait.gene.lm))["EXPR","Pr(>F)"]
    
    # Compute % variance explained by each PC
    trait.var <- var(trait.vals)
    resid.var <- var(trait.gene.lm$residuals)
    trait.gene.propvar <- max(0,(trait.var - resid.var)/trait.var)
    
    # Return as a vector
    return(c(Perc.Var=trait.gene.propvar, ANOVA.PVal=trait.gene.anova.pval, Coeff=trait.gene.coeff))
  })
  names(expr.reg) <- row.names(expr.subset)
  
  # DEPRECATED - BUT THIS MIGHT BE USEFUL TO INCLUDE?
  # Perform correlation tests (Spearman method)
  # expr.cor <- lapply(row.names(expr.subset), function(x){
  #  cor.test(as.numeric(expr.subset[x,]), trait.vals, method="spearman")
  # })
  # names(expr.cor) <- row.names(expr.subset)
  
  expr.trait.reg <- data.frame(row.names=row.names(expr.table))
  expr.trait.reg[,paste0(trait,".Perc.Var")] <- unlist(lapply(expr.reg, "[", "Perc.Var"))
  expr.trait.reg[,paste0(trait,".ANOVA.PVal")] <- unlist(lapply(expr.reg, "[", "ANOVA.PVal"))
  expr.trait.reg[,paste0(trait,".ANOVA.PVal.FDR")] <- p.adjust(expr.trait.reg[,paste0(trait,".ANOVA.PVal")], method="BH")
  
  cat("Distribution of % Var values:\n")
  cat(min(expr.trait.reg[,1]), boxplot.stats(expr.trait.reg[,1])$stats, max(expr.trait.reg[,1]), "\n")
  cat(sum(expr.trait.reg[,3] <= my.args$FDR),"expr features have corrected ANOVA p-value <=", my.args$FDR, "\n")
  
  trait.outdir <- paste0(my.args$OUTDIR, trait, "/")
  if(!file.exists(trait.outdir)) {
    dir.create(path=trait.outdir, recursive = T, showWarnings=F)
  }
  
  trait.outfile <- paste0(trait.outdir, my.args$OUTPUT, "_traitreg.txt")
  
  # Write out to file
  cat("Writing table of all results to:", trait.outfile, "\n")
  if(file.exists(trait.outfile)) {
    cat("WARNING: Overwriting existing file!\n")
  }
  expr.trait.reg <- cbind(GENE=row.names(expr.trait.reg), expr.trait.reg)
  write.table(expr.trait.reg, trait.outfile, sep="\t", row.names=F, quote=F)
  cat("\n")
}
cat("\n")

# Summarize [DEPRECATED - Now done per trait]
# cat("Distribution of % Var values:\n")
# all.perc.var <- c(expr.trait.reg[,grep("[.]Perc[.]Var$", colnames(expr.trait.reg), value=T)], recursive=T)
# cat(min(all.perc.var), boxplot.stats(all.perc.var)$stats, max(all.perc.var), "\n")
# all.pval.FDR <- c(expr.trait.reg[,grep("[.]ANOVA[.]PVal[.]FDR$", colnames(expr.trait.reg), value=T)], recursive=T)
# cat(sum(all.pval.FDR <= my.args$FDR),"expr,trait pairs have corrected ANOVA p-value <=", my.args$FDR, "\n")
# 
# # Write out to file
# cat("Writing table of all results to:", my.args$CORRFILE, "\n")
# if(!file.exists(my.args$OUTDIR)) {
#   dir.create(path=my.args$OUTDIR, recursive = T, showWarnings=F)
# }
# if(file.exists(my.args$CORRFILE)) {
#   cat("WARNING: Overwriting existing file!\n")
# }
# expr.trait.cor <- cbind(GENE=row.names(expr.trait.cor), expr.trait.cor)
# write.table(expr.trait.cor, my.args$CORRFILE, sep="\t", row.names=F, quote=F)

cat("\nScript completed successfully!\n\n")
print(proc.time())
