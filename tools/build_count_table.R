#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 4/18/16
# Adapted from same script in DGRP_Baseline_RNA_Align
#
# Script to build raw count tables
# Rows correspond to Genes/Meta-counts, Columns correspond to samples (aggregated from replicate sequence libraries)
#
# USAGE:
# build_count_table.R [OPTIONS]
# All OPTIONS take the form of PARAM=VALUE (no spaces around '='), the PARAMs are:
#  STUB=      Sets the "stub" part of the count file name, default: STAR
#  OUTPUT=    Sets the output file name, default: combined_samples_gene_counts.txt
#  SAMPLES=   The table containing information about all samples, default: sample_master_table.txt
#  LIBS=      The table containing information about all libraries, default: library_master_table.txt
#  COUNTPATH= The subdirectory containing all HTSeq count output, default: htseq/
#  QCPATH=    The subdirectory to store QC plots, default: sample_QC_figures/
#
# The script looks for SAMPLES and LIBS tables in the current directory (or user must provide full paths)
# The count files are expected to match paths like this:
# [COUNTPATH]/[BATCH]/[LIBNAME]_[STUB]_counts.txt
# (Where BATCH and LIBNAME come from the LIBS table)
# The output file will be saved under [COUNTPATH]/ as well.
#
# Additional QC plots are also output to [QCPATH]/[STUB]_....pdf
#

options(stringsAsFactors=F)

# --- PARSE COMMAND LINE PARAMS --- #

usageStr <- "build_count_table.R [OPTIONS]\n All OPTIONS take the form of PARAM=VALUE (no spaces around '='), the PARAMs are:\n  STUB=\t\t\tSets the \"stub\" part of the count file name, default: STAR\n  OUTPUT=\t\t\tSets the output file name, default: combined_samples_gene_counts.txt\n  SAMPLES=\t\tThe table containing information about all samples, default: sample_master_table.txt  LIBS=\t\t\tThe table containing information about all libraries, default: library_master_table.txt\n  COUNTPATH=\tThe subdirectory containing all HTSeq count output, default: htseq/\n  QCPATH=\t\t\tThe subdirectory to store QC plots, default: sample_QC_figures/\n"
my.args <- commandArgs(trailingOnly=T)
# Debug/Interactive:
# my.args <- c("STUB=STAR_novel", "OUTPUT=combined_samples_novel_counts.txt")

if(length(my.args) > 0) {
  names(my.args) <- unlist(lapply(strsplit(my.args, split="=", fixed=T), "[", 1))
  my.args <- unlist(lapply(strsplit(my.args, split="=", fixed=T), "[", 2))
}

cat("Running with user-specified parameters:\n")
for(param in names(my.args)) {
  cat(param, "=", my.args[param], "\n", sep="")
}

# Sub-routine to set default value
setDefault <- function(param, default) {
  if(!(param %in% names(my.args))) {
    cat("Setting ", param, "=", default, " by default\n", sep="")
    my.args[param] <<- default
  }
}

setDefault("STUB", "STAR")
setDefault("OUTPUT", "combined_samples_gene_counts.txt")
setDefault("SAMPLES", "sample_master_table.txt")
setDefault("LIBS", "library_master_table.txt")
setDefault("COUNTPATH", "htseq/")
setDefault("QCPATH", "sample_QC_figures/")
cat("\n")

# If COUNTPATH or QCPATH is missing trailing "/" add it on
if((!grepl("/$", my.args["COUNTPATH"])) & (my.args["COUNTPATH"] != "")) {
  my.args["COUNTPATH"] <- paste0(my.args["COUNTPATH"], "/")
}

if((!grepl("/$", my.args["QCPATH"])) & (my.args["QCPATH"] != "")) {
  my.args["QCPATH"] <- paste0(my.args["QCPATH"], "/")
}


# --- Load Input Files --- #

# Load the sample table
cat("Loading sample annotation info:",my.args["SAMPLES"],"\n")
sample.table <- read.table(my.args["SAMPLES"], header=T, sep="\t", row.names=1)
cat("\n")

# Load the library table (has batch data, needed for pulling out the HTSeq-count files...)
cat("Loading library info:",my.args["LIBS"],"\n")
library.table <- read.table(my.args["LIBS"], header=T, sep="\t", row.names=1)

# Loop over each sample, load counts from each library
# When there's more than one library, check pairwise correlations
sample.lib.map <- strsplit(sample.table$LIBRARIES, split=",", fixed=T)
names(sample.lib.map) <- row.names(sample.table)
cat("Aggregating read counts from",length(c(sample.lib.map, recursive=T)),"libraries into table for",length(sample.lib.map),"samples.\n")
sample.count.table <- NULL
library.corr.table <- data.frame(SAMPLE=character(0), LIB1=character(0), LIB2=character(0), SCC=numeric(0))
for(smp in row.names(sample.table)) {
  smp.lib.counts <- NULL
  for(lib in sample.lib.map[[smp]]) {
    stopifnot(lib %in% row.names(library.table))
    lib.count.file <- paste0(my.args["COUNTPATH"], library.table[lib,"BATCH"], "/", lib, "_", my.args["STUB"],"_counts.txt")
    lib.counts <- read.table(lib.count.file, header=F, sep="\t", row.names=1)
    if(is.null(smp.lib.counts)) {
      smp.lib.counts <- lib.counts
    } else {
      smp.lib.counts <- cbind(smp.lib.counts, lib.counts)
    }
  }
  colnames(smp.lib.counts) <- sample.lib.map[[smp]]
  # If more than 1 column, check all parwise correlations (spearman)
  if(ncol(smp.lib.counts) > 1) {
    for(i in 1:(ncol(smp.lib.counts)-1)) {
      for(j in (i+1):ncol(smp.lib.counts)) {
        lct.row <- nrow(library.corr.table)+1
        library.corr.table[lct.row,"SAMPLE"] <- smp
        library.corr.table[lct.row,"LIB1"] <- colnames(smp.lib.counts)[i]
        library.corr.table[lct.row,"LIB2"] <- colnames(smp.lib.counts)[j]
        library.corr.table[lct.row,"PCC"] <- cor.test(smp.lib.counts[,i], smp.lib.counts[,j], method="pearson")$estimate
        library.corr.table[lct.row,"SCC"] <- cor.test(smp.lib.counts[,i], smp.lib.counts[,j], method="spearman")$estimate
      }
    }
  }
  # Then sum across columns, add this new column to sample.count.table
  smp.counts <- apply(smp.lib.counts, 1, sum)
  sample.count.table <- cbind(sample.count.table, smp.counts)
}
colnames(sample.count.table) <- row.names(sample.table)

cat("DONE: Tabulated",nrow(sample.count.table),"feature counts for",ncol(sample.count.table),"samples.\n\n")
# EVERYTHING GOES IN SINGLE TABLE, WITH FLAG COLUMN NOW!
# table.file <- "htseq/combined_samples_full_counts.txt"
# cat("Writing to",table.file,"\n")
# write.table(cbind(FEATURE=row.names(sample.count.table), sample.count.table), table.file, row.names=F, sep="\t", quote=F)
               
# Check the library-wise correlations - as long as these are all very high, there is no need to look into library-level batch effects
cat("Checking correlations for re-sequenced libraries:\n\n")
cat("Range of Pearson Correlation Coefficients:\n")
print(range(library.corr.table$PCC))
cat("Boxplot stats of Pearson Correlation Coefficients:\n")
print(boxplot.stats(library.corr.table$PCC)$stats)
cat("Range of Spearman Correlation Coefficients:\n")
print(range(library.corr.table$SCC))
cat("Boxplot stats of Spearman Correlation Coefficients:\n")
print(boxplot.stats(library.corr.table$SCC)$stats)
cat("\n")

# --- U R HERE - 5/18/16 --- #
# NEED A PARAM FOR THESE PLOT FILE...

# If QCPATH doesn't exist, create it:
if((my.args["QCPATH"] != "") & !file.exists(my.args["QCPATH"])) {
  dir.create(my.args["QCPATH"], recursive=T)
}

# Generate bar plots of % reads in __no_feature and __ambiguous categories
cat("Checking % of uniquely aligned reads not assigned to gene features.\n")
perc.ambiguous <- sample.count.table["__ambiguous",] *100 / apply(sample.count.table, 2, sum)
perc.ambiguous <- sort(perc.ambiguous, decreasing=F)
cat("Range of % ambiguously assigned reads across samples:", range(perc.ambiguous),"\n")
cat("Boxplot stats:", boxplot.stats(perc.ambiguous)$stats, "\n")
plot.file <- paste0(my.args["QCPATH"], my.args["STUB"], "_sample_ambiguous_perc.pdf")
cat("Drawing",plot.file,"\n")
pdf(plot.file)
  barplot(perc.ambiguous, border=NA, col="grey", names=NA, ylab="% of Uniquely Aligned Reads w/ Ambiguous Assignment", xlab="Samples ordered by % Ambiguous", main="Ambiguous Read Assignment (HTSeq-count)")
dev.off()

perc.nofeature <- sample.count.table["__no_feature",] *100 / apply(sample.count.table, 2, sum)
perc.nofeature <- sort(perc.nofeature, decreasing=F)
cat("Range of % unassigned reads across samples:", range(perc.nofeature),"\n")
cat("Boxplot stats:", boxplot.stats(perc.nofeature)$stats, "\n")
plot.file <- paste0(my.args["QCPATH"], my.args["STUB"], "_sample_nofeature_perc.pdf")
cat("Drawing",plot.file,"\n")
pdf(plot.file)
  barplot(perc.nofeature, border=NA, col="grey", names=NA, ylab="% of Uniquely Aligned Reads w/ No Feature Assignment", xlab="Samples ordered by % Unassigned", main="Unassigned Reads (HTSeq-count)")
dev.off()

# Now drop all annotation values
# SKIP: Just flag them instead
# drop.rows <- grep("^__", row.names(sample.count.table), value=T)
# cat("Dropping meta-count features:",drop.rows,"\n")
# sample.count.table <- sample.count.table[setdiff(row.names(sample.count.table), drop.rows),]
# cat("\n")

# Flag rRNA features and anything on Mt
ribo.genes <- intersect(row.names(sample.count.table), read.table("~/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-rRNA-genes.txt")[,1])
mito.genes <- intersect(row.names(sample.count.table), read.table("~/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-mt-genes.txt")[,1])
if(length(c(ribo.genes, mito.genes)) > 0) {
  cat("Checking % of reads in ribosomal and mitochondrial features.\n")

  perc.ribo <- apply(sample.count.table[ribo.genes,], 2, sum) * 100 / apply(sample.count.table, 2, sum)
  perc.mito <- apply(sample.count.table[mito.genes,], 2, sum) * 100 / apply(sample.count.table, 2, sum)

  perc.ribo <- sort(perc.ribo, decreasing=F)
  cat("Range of % ribosomal reads across samples:", range(perc.ribo),"\n")
  cat("Boxplot stats:", boxplot.stats(perc.ribo)$stats, "\n")
  plot.file <- paste0(my.args["QCPATH"], my.args["STUB"], "_sample_ribo_perc.pdf")
  cat("Drawing",plot.file,"\n")
  pdf(plot.file)
    barplot(perc.ribo, border=NA, col="grey", names=NA, ylab="% of Reads Assigned to rRNA", xlab="Samples ordered by % Ribosomal", main="Ribosomal Reads")
  dev.off()
  # This is extremely low, nothing to worry about here

  perc.mito <- sort(perc.mito, decreasing=F)
  cat("Range of % mitochondrial reads across samples:", range(perc.mito),"\n")
  cat("Boxplot stats:", boxplot.stats(perc.mito)$stats, "\n")
  plot.file <- paste0(my.args["QCPATH"], my.args["STUB"], "_sample_mito_perc.pdf")
  cat("Drawing",plot.file,"\n")
    pdf(plot.file)
    barplot(perc.mito, border=NA, col="grey", names=NA, ylab="% of Reads Assigned to Mt genes", xlab="Samples ordered by % Mitochondrial", main="Mitochondrial Reads")
  dev.off()
  # A few samples have much higher (~10x) Mt content, but these are still ~1% of total reads, so I don't think there's anything to be concerned about here
}

# Construct flag column:
# Flag mito genes first, so that ribo flag takes precedent
flag.col <- rep("OK", length=nrow(sample.count.table))
names(flag.col) <- row.names(sample.count.table)
cat("Flagging",length(ribo.genes),"rRNA and",length(setdiff(mito.genes, ribo.genes)),"Mt features.\n\n")
flag.col[mito.genes] <- "MT"
flag.col[ribo.genes] <- "RIBO"
# Flag meta counts as META
flag.col[grep("^__", names(flag.col))] <- "META"

# TO DO: Could flag some basic low/rare expression features here, but currently handling that in the post-processing scripts

table.file <- paste0(my.args["COUNTPATH"], my.args["OUTPUT"])
cat("Writing counts for",nrow(sample.count.table),"features across",ncol(sample.count.table),"samples to",table.file,"\n")
write.table(cbind(GENE=row.names(sample.count.table), FLAG=flag.col, sample.count.table), table.file, row.names=F, sep="\t", quote=F)

cat("\nScript completed successfully!\n\n")
print(proc.time())
