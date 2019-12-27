#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 3/8/17
#
# eQTL_tabulate.R
#
# Goal: Tabulate plink results for primary and permuted expression profiles
# 
# Usage: Rscript eQTL_tabulate.R TARGETDIR
#
# TARGETDIR should be the top-level directory for plink results.
# This directory should contain a FEATUREID.qassoc file for each feature
# There should be a Perm{X} subdirectory with matching *.qassoc files for each feature
# The script will auto-detect the set of mapped features and the number of permutations
# The script will skip over directories with no .qassoc files (in case new Perm directories were added later)
#

options(stringsAsFactors=F)

# --- Command Line Param --- #
usageStr="USAGE:\nRscript eQTL_tabulate.R TARGETDIR"
my.args <- commandArgs(trailingOnly=T)
# Interactive testing:
# my.args <- c("microbiome_species_filtered/plink/")

if(length(my.args) != 1) {
  stop("This script takes exactly one parameter.", usageStr)
}

target.dir <- my.args[1]

# If target.dir is missing "/" at the end, add it
if(!grepl("[/]$", target.dir)) target.dir <- paste0(target.dir, "/")

# Make sure target.dir exists
if(!file.exists(target.dir)) {
  stop(target.dir, " does not exist!")
}

cat("Summarizing results in", target.dir, "\n\n")

# Get the list of Perm directories
perm.dirs <- dir(path=target.dir, pattern="^Perm[0-9]+$", include.dirs=T)
cat("Detected", length(perm.dirs), "permutations.\n")
perm.ids <- as.integer(sub("^Perm","",perm.dirs))
perm.dirs <- perm.dirs[order(perm.ids, decreasing=F)]
perm.ids <- sort(perm.ids, decreasing=F)
if(length(perm.ids) == 0) {
  cat("No permutations detected, processing single directory only.\n")
} else if(length(perm.ids) != max(perm.ids)) {
  cat("WARNING: Permutation numbers are not consecutive, missing:\n")
  cat(setdiff(paste0("Perm",1:max(perm.ids)), perm.dirs), "\n")
}

cat("\n")
proc.dirs <- c(target.dir, paste0(rep(target.dir, times=length(perm.dirs)), perm.dirs, "/"))
for(my.dir in proc.dirs) {
  
  # Get the list of features based on available .qassoc files
  qassoc.files <- dir(path=my.dir, pattern=".*[.]qassoc$")
  
  if(length(qassoc.files)==0) {
    cat("Skipping", my.dir, " - no .qassoc files, may have been processed already?\n\n")
  } else {
    # Split based on prefix - each prefix is a different model (e.g. F, M, SexAvg)
    models <- unique(unlist(lapply(strsplit(qassoc.files, split=".", fixed=T), "[", 1)))
    model.qassoc.files <- lapply(models, function(x){grep(paste0("^",x,"[.]"), qassoc.files, value=T)})
    names(model.qassoc.files) <- models
    for(mod in models) {
      names(model.qassoc.files[[mod]]) <- sub(paste0("^",mod,"[.]"),"",sub("[.]qassoc$","",model.qassoc.files[[mod]]))
    }
    cat("Processing", length(models), "models in", my.dir, "\n")

    # Main Processing Loop
    for(mod in models) {
      cat(" Processing model", mod, "with", length(model.qassoc.files[[mod]]), "features.\n")
      
      # Loop over each feature
      # TODO: This loop can be easily parallelized!
      mod.table <- data.frame(FEATURE=character(0), SNPS=character(0), PVALS=character(0))
      mod.row <- 1
      for(feature in names(model.qassoc.files[[mod]])) {
        # TEMP debug code:
        # cat("Processing", feature, "... ")
        
        # Load the qassoc file for this feature
        obs.pval.file <- paste0(my.dir, model.qassoc.files[[mod]][feature])
        if(!file.exists(obs.pval.file)) {
          stop("Missing ", obs.pval.file)
        }
        obs.pval <- read.table(obs.pval.file, header = TRUE, as.is = TRUE)
        
        mod.table[mod.row,"FEATURE"] <- feature
        mod.table[mod.row,"SNPS"] <- ""
        mod.table[mod.row,"PVALS"] <- ""
        
        if (nrow(obs.pval) > 0) {
          # Sort by p-value with most signif on top, 
          # Keep only the SNP ID and P-value columns
          obs.pval <- obs.pval[order(obs.pval[, "P"]), c("SNP", "P")]
          
          mod.table[mod.row,"SNPS"] <- paste(obs.pval$SNP, collapse=",")
          mod.table[mod.row,"PVALS"] <- paste(obs.pval$P, collapse=",")
        }
        
        mod.row <- mod.row + 1
      }
      
      out.file <- paste0(my.dir, mod, ".qassoc.results.txt")
      cat("Writing results to:", out.file, "\n")
      if(file.exists(out.file)) {
        cat("WARNING:", out.file, "already exists - Overwriting!\n")
      }
      
      write.table(mod.table, out.file, sep = "\t", quote = FALSE, row.names = FALSE)
      
      cat("\n")
    }
  }
  
  cat("Now safe to remove all *.qassoc files in", my.dir, "\n\n")
}

cat("Script completed successfully!\n\n")
print(proc.time())
