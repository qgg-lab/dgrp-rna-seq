#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 1/31/18
#
# check_eQTL_pvals.R
#
# Goal: Just examine the p-value distribution for eQTLs in a particular FDR-filtered 
#       result table, make sure they do not come too close to the p-value cutoff used 
#       in the underlying plink analysis.
# 
# Usage: Rscript check_eQTL_pvals.R eqtl_FDR_table.txt
#

options(stringsAsFactors=F)

# --- Command Line Param --- #
usageStr="USAGE:\nRscript check_eQTL_pvals.R eqtl_FDR_table.txt"
my.args <- commandArgs(trailingOnly=T)
# Interactive testing:
# my.args <- c("known_all_novel_genes/plink/ExprF.0.05.fdr.results.txt")

# Fixing the p-value threshold to 1E-6 for now
# TO DO: Make this a second optional parameter
plink.cutoff <- 10 ** -5

if(length(my.args) != 1) {
  stop("This script requires exactly one parameter. ", usageStr)
} else {
  eqtl.file <- my.args[1]
}

# Load the eQTL file
cat("Loading eQTL results from:", eqtl.file, "\n")
eqtl.table <- read.table(eqtl.file, sep="\t", header=T)
cat("Loaded eQTL results for", nrow(eqtl.table), "features.\n\n")

# Now parse out the individual p-values
cat("Extracting individual SNP p-values...\n")
all.pvals <- strsplit(eqtl.table$SIGPVALS, split=",", fixed=T)
all.pvals <- unlist(all.pvals)
all.pvals <- all.pvals[!is.na(all.pvals)]
all.pvals <- as.numeric(all.pvals)
stopifnot(sum(is.na(all.pvals))==0)
cat("Extracted p-values for", length(all.pvals), "individual SNPs.\n")
cat("P-values range from", min(all.pvals), "to", max(all.pvals), "\n")

# Make sure no SNPs are > plink.cutoff
if(any(all.pvals > plink.cutoff)) {
  cat("WARNING:", sum(all.pvals > plink.cutoff), "SNP p-values are above default plink cutoff, may have been changed?")
}

# How many are within 1 order of magnitude of the cutoff?
proximal.cutoff <- plink.cutoff / 10
cat(sum(all.pvals >= proximal.cutoff), "SNP p-values are within 1 order of magnitude of plink cutoff.\n")
cat(" =", mean(all.pvals >= proximal.cutoff)*100, "% of SNPs.\n")
