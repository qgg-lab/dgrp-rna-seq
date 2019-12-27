#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 8/31/16
# Adapted from build_count_table
#
# Script to build raw count tables for individual libraries (NOT combined by sample)
# Rows correspond to Genes/Meta-counts, Columns correspond to libraries
#
# USAGE:
# build_count_table.R [OPTIONS]
# All OPTIONS take the form of PARAM=VALUE (no spaces around '='), the PARAMs are:
#  STUB=      Sets the "stub" part of the count file name, default: STAR
#  OUTPUT=    Sets the output file name, default: library_gene_counts.txt
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
# Additional QC plots are also output to [QCPATH]/library_[STUB]_....pdf
#

options(stringsAsFactors=F)

# --- PARSE COMMAND LINE PARAMS --- #

usageStr <- "build_lib_count_table.R [OPTIONS]\n All OPTIONS take the form of PARAM=VALUE (no spaces around '='), the PARAMs are:\n  STUB=\t\t\tSets the \"stub\" part of the count file name, default: STAR\n  OUTPUT=\t\t\tSets the output file name, default: library_gene_counts.txt\n  LIBS=\t\t\tThe table containing information about all libraries, default: library_master_table.txt\n  COUNTPATH=\tThe subdirectory containing all HTSeq count output, default: htseq/\n  QCPATH=\t\t\tThe subdirectory to store QC plots, default: sample_QC_figures/\n"
my.args <- commandArgs(trailingOnly=T)
# Debug/Interactive:
# my.args <- c("STUB=STAR", "OUTPUT=library_gene_counts.txt")
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
setDefault("OUTPUT", "library_gene_counts.txt")
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

# Load the library table (has batch data, needed for pulling out the HTSeq-count files...)
cat("Loading library info:",my.args["LIBS"],"\n")
library.table <- read.table(my.args["LIBS"], header=T, sep="\t", row.names=1)

# Loop over each library, load counts
cat("Tabulating read counts from",nrow(library.table),"libraries into table format.\n")
library.count.table <- NULL
for(lib in row.names(library.table)) {
  lib.count.file <- paste0(my.args["COUNTPATH"], library.table[lib,"BATCH"], "/", lib, "_", my.args["STUB"],"_counts.txt")
  lib.counts <- read.table(lib.count.file, header=F, sep="\t", row.names=1)
  colnames(lib.counts) <- lib
  if(is.null(library.count.table)) {
    library.count.table <- lib.counts
  } else {
    library.count.table <- cbind(library.count.table, lib.counts)
  }
}

cat("DONE: Tabulated",nrow(library.count.table),"feature counts for",ncol(library.count.table),"libraries.\n\n")

# U R HERE - 8/31/16
# Can probably skip FLAG column for now and just output the table for the other QC stuff I want to do

# Flag rRNA features and anything on Mt
ribo.genes <- intersect(row.names(library.count.table), read.table("~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-rRNA-genes.txt")[,1])
mito.genes <- intersect(row.names(library.count.table), read.table("~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-mt-genes.txt")[,1])

# Construct flag column:
# Flag mito genes first, so that ribo flag takes precedent
flag.col <- rep("OK", length=nrow(library.count.table))
names(flag.col) <- row.names(library.count.table)
cat("Flagging",length(ribo.genes),"rRNA and",length(setdiff(mito.genes, ribo.genes)),"Mt features.\n\n")
flag.col[mito.genes] <- "MT"
flag.col[ribo.genes] <- "RIBO"
# Flag meta counts as META
flag.col[grep("^__", names(flag.col))] <- "META"

# TO DO: Could flag some basic low/rare expression features here, but currently handling that in the post-processing scripts

table.file <- paste0(my.args["COUNTPATH"], my.args["OUTPUT"])
cat("Writing counts for",nrow(library.count.table),"features across",ncol(library.count.table),"libraries to",table.file,"\n")
write.table(cbind(GENE=row.names(library.count.table), FLAG=flag.col, library.count.table), table.file, row.names=F, sep="\t", quote=F)

cat("\nScript completed successfully!\n\n")
print(proc.time())
