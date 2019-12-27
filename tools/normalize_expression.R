#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 10/15/15
#
# USAGE:
# Rscript normalize_expression.R COUNTS=my_counts.txt [OUTDIR=subdir/ SIZE=LIBSZ GENES=gene-info.txt SAMPLES=sample_master_table.txt]
#  COUNTS=  specifies the main count file to normalize
#  OUTDIR=  specifies where to put normalized files (defaults to same path as COUNTS)
#  SIZE=    column in sample table with read counts to normalize to (if not specified, use column totals)
#  GENES=   path to file of gene info, must contain LENGTH column for FPKM normalization (if row IDs don't match up completely, skip FPKM normalization)
#  SAMPLES= table of sample info (defaults to sample_master_table.txt)
#  NORM=    specifies an EdgeR normalization method to apply to library sizes at the RPM step (now defaults to "TMM")
#  CUTOFF=  a minimum normalized expression value, flag genes that are below cutoff in all samples 
#           defaults to NONE for no cutoff - leave flags as they are
#           set to FIT to use an EM model to fit an appropriate cutoff
#  GROUP=  Which column to use for grouping for both EM model and EdgeR normalization
#			Defaults to GROUP, but can specify another column here
#			If the column is missing it will be filled in with LINE_SEX
#
# NOTE: Both the cutoff and EdgeR normalization are applied per "group" 
# This can be specified by adding a GROUP column to sample_master_table.txt
# If no GROUP column is present, then it is filled in with LINE_SEX
#
# NOTE: If counts.txt contains a FLAG column, this will be used to remove genes not marked OK/LOW/RARE before other normalization
#
# This script performs several strategies to produce multiple normalized versions of the output file
# 1) Simple RPM normalization (using column totals)
# 2) Simple FPKM normalization (Normalizing RPM to gene/feature length)
# 3) Simple RLE normalization (Log2 of FPKM+pseudo, subtract feature-wise mean)
#
# TO DO: Round values before outputing the normalized files, otherwise they are much larger than the original count file!
# TO DO: Could put some QC figure generation here?
#
# Output files will have similar name with _counts.txt changed to _rpm.txt, _fpkm.txt, and _rle.txt
# Will NOT overwrite existing files, delete existing before re-running
# If the input file does not end in _counts.txt, then the new suffixes just get appended to the entire input file name
# Output files maintain the 
#


# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

library(edgeR)
library(mixtools)
options(stringsAsFactors=F)

# -- Fixed Parameters -- #
# These could be made adjustable in the command-line args
#
# Change pseudo-count here, this sets the minimum FPKM values when converting to RLE
# TO DO: Automate adapting this value to other data sets, or make it user adjustable
pseudo.count <- 0.001

# Seed random number generation to make the EM mixture model step stable
# TO DO: Allow user adjustment of this:
my.seed <- as.integer(8210678)
cat("Seeding RNG with", my.seed, "\n")
set.seed(my.seed)


# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript normalize_expression.R COUNTS=my_counts.txt [OUTDIR=subdir/ SIZE=LIBSZ GENES=gene-info.txt SAMPLES=sample_master_table.txt CUTOFF=FIT]"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("COUNTS=known_all_novel_genes/combined_samples_known_novel_counts.txt", "OUTDIR=known_all_novel_genes", "GENES=known_all_novel_genes/combined-gene-info.txt", "CUTOFF=FIT")
if(length(my.args) == 0) {
  stop("Requires at least one argument. ", usageStr)
} else {
  cat("Running with params:",my.args,"\n")
}
# Parse param names and values
argSplit=strsplit(my.args, split="=", fixed=T)
my.args <- unlist(lapply(argSplit, "[", 2))
names(my.args) <- unlist(lapply(argSplit, "[", 1))

# Make sure COUNTS param exists
if(!("COUNTS" %in% names(my.args))) {
  stop("Missing COUNTS parameter. ", usageStr)
}

# Fill in default value for missing params
setDefault <- function(param, default) {
  if(!(param %in% names(my.args))) {
    cat("Setting ",param,"=",default," by default\n", sep="")
    my.args[param] <<- default
  }
}

# Make sure COUNTS file exists
if(!file.exists(my.args["COUNTS"])) {
  stop("Missing COUNTS file: ", my.args["COUNTS"])
}

setDefault("OUTDIR", paste0(dirname(my.args["COUNTS"]), "/"))
setDefault("SIZE", "")
setDefault("GENES", "~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-transcriptome-r5.57-gene-info.txt")
setDefault("SAMPLES", "sample_master_table.txt")
setDefault("NORM", "TMM")
setDefault("CUTOFF", "NONE")
setDefault("GROUP", "GROUP")
cat("\n")

# Make sure OUTDIR ends with a "/", set to same directory as COUNTS file if blank
if(my.args["OUTDIR"] == "") {
  my.args["OUTDIR"] <- paste0(dirname(my.args["COUNTS"]), "/")
}
if(!grepl("[/]$", my.args["OUTDIR"])) {
  my.args["OUTDIR"] <- paste0(my.args["OUTDIR"], "/")
}
if(!file.exists(my.args["OUTDIR"])) {
  dir.create(my.args["OUTDIR"], showWarnings=F, recursive=T)
}

# Make sure GENES and SAMPLES files exist (GENES is not required, so just warn if missing)
if(!file.exists(my.args["GENES"])) {
  cat("WARNING:", my.args["GENES"], "does not exist, will skip FPKM normalization\n")
  my.args["GENES"] <- ""
}

if(!file.exists(my.args["SAMPLES"])) {
  stop("Missing SAMPLES file: ", my.args["SAMPLES"])
}

# If CUTOFF is NONE or FIT, leave as character, otherwise try to convert to numeric
# If conversion produces NA, default to NONE with warning
if(!(my.args["CUTOFF"] %in% c("NONE","FIT"))) {
  if(is.na(as.numeric(my.args["CUTOFF"]))) {
    cat("WARNING: Did not recognize CUTOFF parameter, defaulting to NONE.\n")
    my.args["CUTOFF"] <- "NONE"
  }
}

# Base name of COUNTS file comes from shaving off trailing "_counts.txt"
# (done piece-meal in case it ends in .txt but missing "_" or "counts"...)
my.args["BASENAME"] <- sub("[.]txt$", "", basename(my.args["COUNTS"]))
my.args["BASENAME"] <- sub("counts$", "", my.args["BASENAME"])
my.args["BASENAME"] <- sub("[_]$", "", my.args["BASENAME"])

# Specify output files based on OUTDIR and BASENAME
# DO NOT OVERWRITE ANY OF THESE - Quit with error if they exist already
my.args["RPMFILE"] <- paste0(my.args["OUTDIR"], my.args["BASENAME"], "_rpm.txt")
if(file.exists(my.args["RPMFILE"])) {
  stop(my.args["RPMFILE"], " already exists - delete all normalized tables before re-running this script!")
}
my.args["FPKMFILE"] <- paste0(my.args["OUTDIR"], my.args["BASENAME"], "_fpkm.txt")
if(file.exists(my.args["FPKMFILE"])) {
  stop(my.args["FPKMFILE"], " already exists - delete all normalized tables before re-running this script!")
}
my.args["RLEFILE"] <- paste0(my.args["OUTDIR"], my.args["BASENAME"], "_rle.txt")
if(file.exists(my.args["RLEFILE"])) {
  stop(my.args["RLEFILE"], " already exists - delete all normalized tables before re-running this script!")
}


# --- LOAD INPUT --- #

# Load the count table
count.table <- read.table(my.args["COUNTS"], header=T, sep="\t", row.names=1)
# If first column == FLAG, use that to filter the table, then drop the column
dropped.features <- 0
if(colnames(count.table)[1] == "FLAG") {
  dropped.features <- sum(!(count.table$FLAG %in% c("OK","RARE","LOW")))
  count.table <- count.table[count.table$FLAG %in% c("OK","RARE","LOW"),]
  gene.flags <- count.table$FLAG
  count.table <- count.table[,2:ncol(count.table)]
} else {
  gene.flags <- rep("OK", times=nrow(count.table))
}
names(gene.flags) <- row.names(count.table)

# Remove 'X' that R automatically puts in front of column names
colnames(count.table) <- sub("^X", "", colnames(count.table))
features <- row.names(count.table)
cat("Loaded counts for",length(features),"features across",ncol(count.table),"samples from",my.args["COUNTS"],"\n")
if(dropped.features > 0) {
  cat("(Dropped",dropped.features,"features based on FLAG column)\n")
}

# Load the master sample table
sample.table <- read.table(my.args["SAMPLES"], header=T, sep="\t", row.names=1)
cat("Loaded info on",nrow(sample.table),"samples from",my.args["SAMPLES"],"\n")

# If SIZE is specified, make sure the column is present in sample.table
if(my.args["SIZE"] != "") {
  if(!(my.args["SIZE"] %in% colnames(sample.table))) {
    cat("WARNING:", my.args["SAMPLES"], "does not contain a", my.args["SIZE"], "column for desired RPM normalization - will use sample count totals instead.\n")
    my.args["SIZE"] == ""
  }
}

# How many samples are in common between the two?
all.samples <- intersect(row.names(sample.table),colnames(count.table))
cat(length(all.samples),"samples matched up between sample info and count table.\n")
if(length(all.samples) < nrow(sample.table)) {
  cat("WARNING: Dropping",length(setdiff(row.names(sample.table),all.samples)),"rows from sample info table.\n")
}
if(length(all.samples) < ncol(count.table)) {
  cat("WARNING: Dropping",length(setdiff(colnames(count.table),all.samples)),"columns from count table.\n")
}
sample.table <- sample.table[all.samples,]
count.table <- count.table[,all.samples]

# Check for GROUP column, default to LINE_SEX if not present
# Either way, let the user know what's happening
# If there's a GROUP column, use that, otherwise just combine LINE_SEX
if(my.args["GROUP"] %in% colnames(sample.table)) {
  cat("Using",my.args["GROUP"],"column to define", length(unique(sample.table[,my.args["GROUP"]])), "sample groups.\n")
} else {
  sample.table[,my.args["GROUP"]] <- paste(sample.table$LINE, sample.table$SEX, sep="_")
  cat("Using LINE_SEX to define", length(unique(sample.table[,my.args["GROUP"]])), "sample groups.\n")
}

# Get feature lengths from file specified by GENES parameter
feature.length <- c()
if(my.args["GENES"] == "") {
  cat("No feature info file specified to GENES parameter - FPKM normalization will be SKIPPED.\n")
} else if(!file.exists(my.args["GENES"])) {
  cat("WARNING: Cannot find", my.args["GENES"], "- FPKM normalization will be SKIPPED.\n")
} else {
  feature.info.table <- read.table(my.args["GENES"], header=T, sep="\t", row.names=1)
  if(!all(features %in% row.names(feature.info.table))) {
    cat("WARNING:", my.args["GENES"], "did not contain info on all count features in", my.args["COUNTS"], "- FPKM normalization will be SKIPPED.\n")
  } else {
    feature.length <- feature.info.table[features,"LENGTH"]
    names(feature.length) <- features
    cat("Loaded feature length for",length(feature.length),"genes/features from",my.args["GENES"],"\n")
  }
}
cat("\n")


# --- RPM Normalization --- #

sample.sizes <- c()
if(my.args["SIZE"] == "") {
  cat("Performing RPM normalization using sample count sums.\n")
  sample.sizes <- apply(count.table, 2, sum)
} else {
  cat("Performing RPM normalization using",my.args["SIZE"],"column as library size.\n")
  sample.sizes <- sample.table[,my.args["SIZE"]]
  names(sample.sizes) <- row.names(sample.table)
}
stopifnot(all(names(sample.sizes) == colnames(count.table)))
  
cat("Sample lib sizes range from", round(min(sample.sizes)/(10**6), digits=1), "M to", round(max(sample.sizes)/(10**6), digits=1), "M\n")
cat("Median sample lib size is", round(median(sample.sizes)/(10**6), digits=1), "M\n")
cat("Sample lib size IQR is", paste(round(boxplot.stats(sample.sizes)$stats[c(2,4)]/(10**6), digits=1), collapse=" - "), "M\n")

# Try using EdgeR normalization here...
if(my.args["NORM"] == "none") {
  rpm.table <- count.table
  for(j in 1:ncol(rpm.table)) {
    rpm.table[,j] <- count.table[,j] * (10**6) / sample.sizes[j]
  }
  cat("\n")
} else {
  dgelist.general <- DGEList(counts=count.table, lib.size=sample.sizes, group=sample.table[,my.args["GROUP"]], genes=feature.info.table[row.names(count.table),])
  cat("Normalizing by", my.args["NORM"], "method built into EdgeR.\n")
  dgelist.general <- calcNormFactors(dgelist.general, method=my.args["NORM"])
  # range(dgelist.general$samples[grep("F", row.names(sample.table), value=T),"norm.factors"])
  # range(dgelist.general$samples[grep("M", row.names(sample.table), value=T),"norm.factors"])
  rpm.table <- as.data.frame(cpm(dgelist.general, normalized.lib.sizes=T))
}

cat("Global median RPM expression value:", round(median(as.matrix(rpm.table)), digits=1), "\n")
cat("Global range of RPM values:", round(min(as.matrix(rpm.table)), digits=1), "to", round(max(as.matrix(rpm.table)), digits=1), "\n")
sample.median.rpm <- apply(rpm.table, 2, median)
cat("Sample-wise median RPM values range from:", round(min(sample.median.rpm), digits=1), "to", round(max(sample.median.rpm), digits=1), "\n")
cat("IQR of sample-wise medians:", paste(round(boxplot.stats(sample.median.rpm)$stats[c(2,4)], digits=1), collapse=" - "), "\n")

# TO DO: Could plot RPM distributions here?

# Write file out
cat("Writing RPM table to", my.args["RPMFILE"], "\n")
write.table(cbind(GENE=row.names(rpm.table), FLAG=gene.flags[row.names(rpm.table)], rpm.table), my.args["RPMFILE"], sep="\t", row.names=F, quote=F)
cat("\n")


# --- FPKM Normalization --- #

if(length(feature.length) == 0) {
  cat("No feature length information - skipping FPKM normalization.\n")
  fpkm.table <- rpm.table
} else {
  cat("Performing FPKM normalizing using feature lengths.\n")

  stopifnot(all(names(feature.length) == row.names(rpm.table)))
  fpkm.table <- rpm.table
  for(j in 1:ncol(fpkm.table)) {
    fpkm.table[,j] <- rpm.table[,j] * (10**3) / feature.length
  }
  
  cat("Global median FPKM expression value:", round(median(as.matrix(fpkm.table)), digits=1), "\n")
  cat("Global range of FPKM values:", round(min(as.matrix(fpkm.table)), digits=1), "to", round(max(as.matrix(fpkm.table)), digits=1), "\n")
  sample.median.fpkm <- apply(fpkm.table, 2, median)
  cat("Sample-wise median FPKM values range from:", round(min(sample.median.fpkm), digits=1), "to", round(max(sample.median.fpkm), digits=1), "\n")
  cat("IQR of sample-wise medians:", paste(round(boxplot.stats(sample.median.fpkm)$stats[c(2,4)], digits=1), collapse=" - "), "\n")

  # Write file out - DO NOT OVERWRITE
  cat("Writing FPKM table to", my.args["FPKMFILE"], "\n")
  write.table(cbind(GENE=row.names(fpkm.table), FLAG=gene.flags[row.names(fpkm.table)], fpkm.table), my.args["FPKMFILE"], sep="\t", row.names=F, quote=F)
}
cat("\n")


# --- (Optional) Redefine threshold for LOW/RARE based on distribution of Log2 avg normalized counts --- #

if(my.args["CUTOFF"] != "NONE") {
  cat("Readjusting LOW/RARE flags based on normalized expression values.\n")
  
  unit.type <- if(length(feature.length) == 0){"RPM"}else{"FPKM"}
  
  if(my.args["CUTOFF"] == "FIT") {
    cat("Using mixture model to determine appropriate cutoff.\n")

    # Get all log2 FPKM values, drop the -Inf values
    all.log2.fpkm <- c(log2(fpkm.table), recursive=T)
    all.log2.fpkm <- all.log2.fpkm[all.log2.fpkm > -Inf]

    # Downsample to some max number of expr values, e.g. 1M
    if(length(all.log2.fpkm) > (10**6)) {
      cat("Down-sampling to 1M random expr values.\n")
      all.log2.fpkm <- sample(all.log2.fpkm, size=10**6, replace=F)
    }
    # Should still get a cutoff point around -1.8

    # Fit mixture model of this distribution of values
    cat("Fitting mixture model on Gene-wise Avg Log2 FPKM")
    mixmdl <- normalmixEM(all.log2.fpkm, fast=T)
    # Figure out which is the lower cluster
    low.clust <- which.min(mixmdl$mu)
    # Which values have a 95% change (posterior >= 0.95) of being in cluster 1?
    low.clust.values <- all.log2.fpkm[mixmdl$posterior[,low.clust] >= 0.95]
    # Cutoff point is largest value that is still strongly in the low cluster
    cutoff.point <- max(low.clust.values)
    # Draw a plot
    pdf(paste0(my.args["OUTDIR"], "normalize_expression_mixmdl_density.pdf"))
    plot(mixmdl, which=2, xlab2=paste("Log2", unit.type), main2="Bimodal Mixture Model")
    abline(v=cutoff.point, col="blue", lwd=2)
    dev.off()
    # Convert back from log2 scale
    cutoff.point <- 2**cutoff.point
    cat("Mixture model determined optimal cutoff point =", cutoff.point, unit.type, "\n")
  } else {
    cutoff.point <- as.numeric(my.args["CUTOFF"])
    if(is.na(cutoff.point)) {
      cat("Did not recognized value for CUTOFF =", my.args["CUTOFF"], "- skipping flag adjustments.\n")
    } else {
      cat("Using user-specified cutoff =", cutoff.point, "\n")
    }
  }
  
  # Apply the cutoff
  # Genes with no single line,sex pair with both replicates >= cutoff.point are marked LOW
  # Genes with <5% of line,sex pairs having both replicates >= cutoff.point are marked RARE
  if(!is.na(cutoff.point)) {
    cat("Reassigning LOW/RARE flags based on cutoff point =", cutoff.point, "\n")
    cat("Each gene or feature must be >=", cutoff.point, "in all samples of at least one group.\n")
    new.flags <- rep("OK", times=nrow(fpkm.table))
    names(new.flags) <- row.names(fpkm.table)
	sample.groups <- sample.table[,my.args["GROUP"]]
    names(sample.groups) <- row.names(sample.table)
    stopifnot(all(names(sample.groups)==colnames(fpkm.table)))
    
    group.min.table <- NULL
    for(grp in unique(sample.groups)) {
      grp.samples <- names(sample.groups[sample.groups==grp])
      if(length(grp.samples)==1) {
        group.min.table <- cbind(group.min.table, fpkm.table[,grp.samples])
      } else {
        group.min.table <- cbind(group.min.table, apply(fpkm.table[,grp.samples], 1, min))
      }
    }
    colnames(group.min.table) <- unique(sample.groups)
    
    # Those which have <5% of max.group.min < cutoff.point marked RARE
    new.flags[apply(group.min.table >= cutoff.point, 1, mean) < 0.05] <- "RARE"
    # Those which have all group.min < cutoff.point marked LOW
    max.group.min <- apply(group.min.table, 1, max)
    new.flags[max.group.min < cutoff.point] <- "LOW"
    
    # Report:
    cat("Original flag counts:\n")
    print(table(gene.flags))
    cat("\nTransitions:\n")
    for(old.flg in unique(gene.flags)) {
      for(new.flg in unique(new.flags)) {
        if(any((gene.flags == old.flg) & (new.flags == new.flg))) {
          cat(sum((gene.flags == old.flg) & (new.flags == new.flg)), old.flg, "->", new.flg, "\n")
        }
      }
    }
    cat("\n")
    cat("New flag counts:\n")
    print(table(new.flags))
    cat("\n")
    
    stopifnot(all(names(gene.flags)==names(new.flags)))
    gene.flags <- new.flags
    
    # Now re-write output file, either RPM or FPKM depending on unit.type
    if(unit.type == "RPM") {
      cat("Re-writing RPM table with updated flags to", my.args["RPMFILE"], "\n")
      write.table(cbind(GENE=row.names(rpm.table), FLAG=gene.flags[row.names(rpm.table)], rpm.table), my.args["RPMFILE"], sep="\t", row.names=F, quote=F)
    } else if(unit.type == "FPKM") {
      cat("Re-writing FPKM table with updated flags to", my.args["FPKMFILE"], "\n")
      write.table(cbind(GENE=row.names(fpkm.table), FLAG=gene.flags[row.names(fpkm.table)], fpkm.table), my.args["FPKMFILE"], sep="\t", row.names=F, quote=F)
    } else {
      stop("Unknown unit.type = ", unit.type)
    }
    
  } else {
    cat("No appropriate cutoff point set, skipping!\n")
  }
  cat("\n")
}

# --- RLE Normalization --- #

cat("Performing RLE normalization using mean log2 expression values for each feature.\n")

# General version - just get % of samples w/ non-zero expression
# TO DO: It might be good to have a parameter to change 0 to some other threshold
min.rep.expr.perc <- apply(fpkm.table, 1, function(x){mean(x > 0)})

non.rep.features <- names(min.rep.expr.perc)[min.rep.expr.perc == 0]

cat("Removing", length(non.rep.features), "features with 0 expression in all samples.\n")
keep.features <- setdiff(row.names(fpkm.table), non.rep.features)
fpkm.table <- fpkm.table[keep.features,]

# Compute minimum non-zero value
min.expr <- apply(fpkm.table, 2, function(x){min(x[x > 0])})
min.expr <- min(min.expr)
cat("Minimum non-zero FPKM value (remaining features):", min.expr, "\n")
cat("Using pseudo-count for log2 transformation", pseudo.count, "\n")
if(min.expr < pseudo.count) {
  cat("WARNING: pseudo-count is larger than minimum non-zero expression value.")
}
cat("\n")

rle.table <- log2(fpkm.table+pseudo.count)
stopifnot(sum(rle.table == -Inf)==0)
stopifnot(sum(rle.table == Inf)==0)
stopifnot(sum(is.na(rle.table))==0)

# Subtract feature-wise means
feature.means <- apply(rle.table, 1, mean)
cat("Feature mean log2 values range from", min(feature.means), "to", max(feature.means), "\n")
for(j in 1:ncol(rle.table)) {
  rle.table[,j] <- rle.table[,j] - feature.means
}
feature.mean.rle <- apply(rle.table, 1, mean)
cat("Feature-wise RLE mean values range from", min(feature.mean.rle), "to", max(feature.mean.rle), "\n")
feature.rle.var <- apply(rle.table, 1, var)
cat("Feature-wise RLE variances range from", min(feature.rle.var), "to", max(feature.rle.var), "\n")

cat("Global median RLE expression value:", round(median(as.matrix(rle.table)), digits=1), "\n")
cat("Global range of RLE values:", round(min(as.matrix(rle.table)), digits=1), "to", round(max(as.matrix(rle.table)), digits=1), "\n")
sample.median.rle <- apply(rle.table, 2, median)
cat("Sample-wise median RLE values range from:", round(min(sample.median.rle), digits=1), "to", round(max(sample.median.rle), digits=1), "\n")
cat("IQR of sample-wise medians:", paste(round(boxplot.stats(sample.median.rle)$stats[c(2,4)], digits=1), collapse=" - "), "\n")

# Write file out
cat("Writing RLE table to", my.args["RLEFILE"], "\n")
write.table(cbind(GENE=row.names(rle.table), FLAG=gene.flags[row.names(rle.table)], rle.table), my.args["RLEFILE"], sep="\t", row.names=F, quote=F)

cat("\n")


cat("Script completed successfully!\n")
print(proc.time())
quit("no")
