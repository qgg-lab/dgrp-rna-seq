#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 6/7/16
#
# trait_correlations.R
#
# Goal: Load in expr feature line means and assess correlation with all available line features
# 
# Usage:
# Rscript trait_correlations [OPTIONS] EXPR=subdir/..._line_means.txt
#  EXPR=  The file containing expression values to assess heritability on
#  LOG2=   Take Log2 of line means before doing correlation against traits
#  SEX=   F/M (Default=infer from EXPR file), determines which sex to compare to for line means
#  TRAITS=  Path to table of trait data for appropiate sex, 
#   defaults to ~/Resources/DGRP/Freeze2/line_means_Wol_(female|male).txt
#
# Automatically puts all output in a traitCorr subdirectory in the same parent directory as data file, with similar basenames
#
# TO DO: Try Benjamini-Hochberg FDR instead of Bonferroni correction?
# TO DO: Should perform permutations to compute non-parametric FDR for the correlations?
# TO DO: Create a parameter to make the Wolbachia adjusted version optional
# TO DO: Create optional parameters to specify pseudo-count [PSEUDO=], 
#   as well as some of the other args already specified below
#

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

options(stringsAsFactors=F)

# -- Fixed Parameters -- #
# These could be made adjustable in the command-line args
pseudo.count <- 0.001
fdr.cutoff <- 0.05
# Path to trait line mean tables (when auto-determining)
trait.dir <- "~/Resources/DGRP/Freeze2/"
# Don't need these any more?
# sample.file <- "sample_master_table.txt"
# line.file <- "~/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/line_wolbachia_inversion_info.txt"

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript trait_correlations.R [OPTIONS] EXPR=subdir/..._line_means.txt"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("EXPR=microbiome_species_filtered/genVar/combined_microbe_filtered_rpm_WolAdj_SexAvg_line_means.txt")
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
  if(grepl("[_]F[_]line[_]means.txt$", my.args$EXPR)) {
    my.args$SEX <- "F"
  } else if(grepl("[_]M[_]line[_]means.txt$", my.args$EXPR)) {
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

setDefault("TRAITS","")
# When traits file is not specified, auto-determine
if(my.args$TRAITS=="") {
  my.args$TRAITS <- paste0(trait.dir, "line_means_Wol_")
  if(my.args$SEX=="F") {
    my.args$TRAITS <- paste0(my.args$TRAITS, "fe")
  }
  my.args$TRAITS <- paste0(my.args$TRAITS, "male.txt")
}

setDefault("LOG2",F)
# If "0" or "1" map these to "FALSE" and "TRUE"
if(is.character(my.args$LOG2) & (my.args$LOG2 == "0")){my.args$LOG2 <- F}
if(is.character(my.args$LOG2) & (my.args$LOG2 == "1")){my.args$LOG2 <- T}
# Convert any other type to logical
if(!is.logical(my.args$LOG2)){my.args$LOG2 <- as.logical(my.args$LOG2)}
# If conversion produced NA, just set to F
if(is.na(my.args$LOG2)) {
  cat("WARNING: Did not recognize value of LOG2 parameter, setting to F by default.\n")
  my.args$LOG2 <- F
}
# UPDATE: 6/20/16
# LOG2 transformation is no longer applied here, it is applied directly to the line_means file
# Generate a big warning for now
# TO DO: this option should be completely removed, along with pseudo-count option!
if(my.args$LOG2) {
  cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("!           WARNING           !\n")
  cat("! LOG2 OPTION IS DEPRECATED   !\n")
  cat("! Are you sure the line means !\n")
  cat("! are not already log scale?  !\n")
  cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
}
  
my.args$PSEUDO <- pseudo.count
my.args$FDR <- fdr.cutoff

# Extract path to EXPR input file, use it as the output dir
# (If EXPR data is under a genVar subdirectory, drop that from OUTDIR path)
my.args$OUTDIR <- paste0(sub("genVar$","", dirname(my.args$EXPR)), "/traitCorr/")

# Base name of COUNTS file comes from shaving off the trailing _[FM]_line_means.txt
my.args$BASENAME <- sub("[.]txt$", "", basename(my.args$EXPR))
my.args$BASENAME <- sub("[_]line[_]means$", "", my.args$BASENAME)
my.args$BASENAME <- sub("[_][FM]$", "", my.args$BASENAME)

# OUTSTUB will be OUTDIR/BASENAME(_Log2)
# This is used to determine output file paths
my.args$OUTSTUB <- paste0(my.args$OUTDIR, my.args$BASENAME)
if(my.args$LOG2) {
  my.args$OUTSTUB <- paste0(my.args$OUTSTUB, "_Log2")
}

# Correlation output table will be OUTSTUB_traitCorr.txt
# Just overwrite if they exist - these output files are tracked by git so it's fine
my.args$CORRFILE <- paste0(my.args$OUTSTUB, "_", my.args$SEX, "_traitCorr.txt")


# --- LOAD INPUT --- #

# Load the line mean table
expr.table <- read.table(my.args$EXPR, header=T, sep="\t", row.names=1)
# Remove 'X' that R automatically puts in front of column names
colnames(expr.table) <- sub("^X", "", colnames(expr.table))
features <- row.names(expr.table)
cat("Loaded counts for",length(features),"features across",ncol(expr.table),"lines from",my.args$EXPR,"\n")

trait.table <- read.table(my.args$TRAITS, header=T, sep="\t", row.names=1)
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

if(my.args$LOG2) {
  cat("Converting expression data to Log2 scale with pseudo-count =",my.args$PSEUDO,"\n")
  expr.table <- log2(expr.table+my.args$PSEUDO)
}
cat("\n")



# --- COMPUTE CORRELATIONS --- #

# This table will store trait correlations and p-values for each gene
expr.trait.cor <- data.frame(row.names=row.names(expr.table))

# loop over each trait and run cor.test to get correlation and p-value
# 12/1/16 - USING SPEARMAN COR (DOES NOT ASSUME BIVARIATE NORMAL DISTRIBUTIONS)
for(trait in row.names(trait.table)) {
  # Extract trait values to a vector, drop the ones that are NA
  trait.vals <- as.numeric(trait.table[trait,])
  names(trait.vals) <- colnames(trait.table)
  trait.vals <- trait.vals[!is.na(trait.vals)]
  cat("Computing expression correlations for",trait,"trait data from",length(trait.vals),"lines.\n")
  expr.subset <- expr.table[,names(trait.vals)]
  # Perform correlation tests (Pearson method)
  expr.cor <- lapply(row.names(expr.subset), function(x){
    cor.test(as.numeric(expr.subset[x,]), trait.vals, method="spearman")
  })
  names(expr.cor) <- row.names(expr.subset)
  expr.trait.cor[,paste0(trait,".cor")] <- unlist(lapply(expr.cor, function(x){x$est}))
  expr.trait.cor[,paste0(trait,".pval")] <- unlist(lapply(expr.cor, function(x){x$p.value}))
  expr.trait.cor[,paste0(trait,".pval.FDR")] <- p.adjust(expr.trait.cor[,paste0(trait,".pval")], method="BH")
}
cat("\n")

# Summarize
cat("Distribution of correlation values:\n")
all.cor <- c(expr.trait.cor[,grep("[.]cor$", colnames(expr.trait.cor), value=T)], recursive=T)
cat(min(all.cor), boxplot.stats(all.cor)$stats, max(all.cor), "\n")
all.pval.FDR <- c(expr.trait.cor[,grep("[.]pval[.]FDR$", colnames(expr.trait.cor), value=T)], recursive=T)
cat(sum(all.pval.FDR <= my.args$FDR),"expr,trait pairs have corrected p-value <=", my.args$FDR, "\n")

# Write out to file
cat("Writing table of all results to:", my.args$CORRFILE, "\n")
if(!file.exists(my.args$OUTDIR)) {
  dir.create(path=my.args$OUTDIR, recursive = T, showWarnings=F)
}
if(file.exists(my.args$CORRFILE)) {
  cat("WARNING: Overwriting existing file!\n")
}
expr.trait.cor <- cbind(GENE=row.names(expr.trait.cor), expr.trait.cor)
write.table(expr.trait.cor, my.args$CORRFILE, sep="\t", row.names=F, quote=F)

cat("\nScript completed successfully!\n\n")
print(proc.time())
