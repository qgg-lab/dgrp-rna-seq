#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 12/14/16
#
# eQTL_perm_fdr.R
#
# Goal: Use primary and permuted eQTL results (plink) to determine FDR thresholds on each expression feature
# 
# Usage: Rscript eQTL_perm_fdr.R TARGETDIR [FDR]
#
# TARGETDIR should be the top-level directory for plink results.
# This directory should contain a MODEL.qassoc.results.txt file for the feature set from each model
# There should be a Perm{X} subdirectory with matching MODEL.qassoc.results.txt files
# The script will auto-detect the group IDs and the number of permutations
# FDR should be > 0 and < 1 (default = 0.05)
#
# TO DO: Add doMC/dopar multi-threading
#

options(stringsAsFactors=F)

# --- Command Line Param --- #
usageStr="USAGE:\nRscript eQTL_perm_fdr.R TARGETDIR [FDR]"
my.args <- commandArgs(trailingOnly=T)
# Interactive testing:
# my.args <- c("microbiome_species_filtered/plink/")

# Default FDR = 5%
fdr <- 0.05

if(length(my.args) < 1) {
  stop("This script requires at least one parameter. ", usageStr)
} else {
  target.dir <- my.args[1]
  if(length(my.args) > 1) {
    # Optional second param is FDR
    fdr <- as.numeric(my.args[2])
    if((fdr <= 0) | (fdr >= 1)) {
      stop("Invalid FDR: ", my.args[2], " - must be between 0 and 1. ", usageStr)
    } else {
      cat("Setting FDR =", fdr, "\n")
    }
    if(length(my.args) > 2) {
      stop("Too many parameters specified. ", usageStr)
    }
  } else {
    cat("Default FDR =", fdr, "\n")
  }
}

# If target.dir is missing "/" at the end, add it
if(!grepl("[/]$", target.dir)) target.dir <- paste0(target.dir, "/")

# Make sure target.dir exists
if(!file.exists(target.dir)) {
  stop(target.dir, " does not exist!")
}

cat("Computing", fdr, "FDR cutoffs for all eQTL results in", target.dir, "\n\n")


# Get the list of Perm directories
perm.dirs <- dir(path=target.dir, pattern="^Perm[0-9]+$", include.dirs=T)
cat("Detected", length(perm.dirs), "permutations.\n")
perm.ids <- as.integer(sub("^Perm","",perm.dirs))
perm.dirs <- perm.dirs[order(perm.ids, decreasing=F)]
perm.ids <- sort(perm.ids, decreasing=F)
if(length(perm.ids) != max(perm.ids)) {
  cat("WARNING: Permutation numbers are not consecutive, missing:\n")
  cat(setdiff(paste0("Perm",1:max(perm.ids)), perm.dirs), "\n")
}


# Get the list of feature groups based on available *.qassoc.results.txt
qassoc.files <- dir(path=target.dir, pattern=".*[.]qassoc[.]results[.]txt$")
models <- sub("[.]qassoc[.]results[.]txt$", "", qassoc.files)
names(qassoc.files) <- models
cat("Processing", length(models), "models:", models, "\n\n")


# Main Processing Loop
for(mod in models) {
  
  cat("Processing model ", mod, ", loading ", target.dir, qassoc.files[mod],"\n", sep="")
  
  mod.table <- read.table(paste0(target.dir,qassoc.files[mod]), header=T, sep="\t")
  cat(mod, "contains", nrow(mod.table), "features.\n")
  
  # Split the mod.table into a list of tables, one per feature, with SNP and PVal columns
  mod.snps <- strsplit(mod.table$SNPS, ",")
  names(mod.snps) <- mod.table$FEATURE
  mod.pvals <- lapply(strsplit(mod.table$PVALS, ","), as.numeric)
  names(mod.pvals) <- mod.table$FEATURE
  # TO DO: Parallelize!
  mod.data <- lapply(mod.table$FEATURE, function(feature){
    data.frame(SNP=mod.snps[[feature]], PVAL=mod.pvals[[feature]])
  })
  names(mod.data) <- mod.table$FEATURE
  
  rm(mod.table, mod.snps, mod.pvals)
  
  # Now load each corresponding table from Perm directories, 
  # combine into single tables for each feature, arranged in list structure much like
  perm.data <- lapply(mod.data, function(x){data.frame(SNP=character(0), PVAL=numeric(0))})
  
  for(pdir in perm.dirs) {
    mod.perm.file <- paste0(target.dir, pdir, "/", qassoc.files[mod])
    cat("Loading permuted results from", mod.perm.file, "\n")
    if(!file.exists(mod.perm.file)) {
      stop("Missing ", mod.perm.file)
    }
    perm.table <- read.table(mod.perm.file, header=T, sep="\t")
    if(!all(perm.table$FEATURE == names(mod.data))) {
      stop(mod.perm.file, " does not contain a matching set of features")
    }
    
    # Split into a list of tables, append to perm.data
    perm.snps <- strsplit(perm.table$SNPS, ",")
    names(perm.snps) <- perm.table$FEATURE
    perm.pvals <- lapply(strsplit(perm.table$PVALS, ","), as.numeric)
    names(perm.pvals) <- perm.table$FEATURE
    
    # TO DO: Parallelize
    for(feature in names(mod.data)) {
      perm.data[[feature]] <- rbind(perm.data[[feature]], data.frame(SNP=perm.snps[[feature]], PVAL=perm.pvals[[feature]]))
    }
  }
  cat("\n")
  
  # Set up FDR results table
  fdr.out <- data.frame(row.names=names(mod.data), FEATURE=names(mod.data))
  fdr.out$FDR <- fdr
  fdr.out$SIGSNPS <- as.character(NA)
  fdr.out$SIGPVALS <- as.character(NA)
  # Loop over each feature
  # TODO: This loop can be easily parallelized!
  for(feature in names(mod.data)) {
    # TEMP debug code:
    # cat("Processing", feature, "... ")
    
    if(nrow(mod.data[[feature]]) > 0) {
      # Make sure SNPs sorted by p-value with most signif on top, 
      mod.data[[feature]] <- mod.data[[feature]][order(mod.data[[feature]]$PVAL),]
      
      # Add an FDR column to this table
      mod.data[[feature]][,"FDR"] <- 0
      
      # TO DO: This can sped up with an lapply command?
      for(j in 1:nrow(mod.data[[feature]])) {
        mod.data[[feature]][j,"FDR"] <- sum(perm.data[[feature]]$PVAL <= mod.data[[feature]][j,"PVAL"])/length(perm.dirs)/j
      }
      
      # identify the cutoff point
      thres <- which(cummax(mod.data[[feature]]$FDR) < fdr)
      if(length(thres) > 0) {
        fdr.out[feature,"FDR"] <- mod.data[[feature]][max(thres), 3]
        fdr.out[feature,"SIGSNPS"] <- paste(mod.data[[feature]][1:max(thres), 1], collapse = ",")
        fdr.out[feature,"SIGPVALS"] <- paste(mod.data[[feature]][1:max(thres), 2], collapse = ",")
      }
      cat(length(thres), "QTLs at", fdr, "FDR.\n")
    } else {
      cat("0 SNPs returned by plink.\n")
    }
  }
  
  out.file <- paste0(target.dir, mod, ".", fdr, ".fdr.results.txt")
  cat("Writing results to:", out.file, "\n")
  if(file.exists(out.file)) {
    cat("WARNING:", out.file, "already exists - Overwriting!\n")
  }
  
  write.table(fdr.out, out.file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n\n")
}

cat("Script completed successfully!\n\n")
print(proc.time())
