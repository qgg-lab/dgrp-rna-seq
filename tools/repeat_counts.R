#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 12/16/15
# 
# GOAL: Collect the read counts from BWA alignment against RepBase
# Summarize the counts across all samples and within samples (percentage) and within taxa
# 
# Usage:
# Rscript repeat_counts.R [OPTIONS]
#  SAMPLES= The table with sample information (default: sample_master_table.txt)
#  LIBRARIES= The table with library batch information (default: library_master_table.txt)
#  GROUPS= The table of class groupings for all target sequences in RepBase DB
#    Default: /home/ljeveret/Resources/RepBase/dmel.repbase-20.01.annot.txt
#  BWAPATH= The path to top-level BWA results where the count files live (default: bwa_repeat/)
#  COUNTS= The name of count files in each BWA output directory (default: target_counts.txt)
#  PRE= The file containing fastq read counts BEFORE repeat filter (default: microbefiltered/fastq.gz.stats)
#  POST= The file containing fastq read counts AFTER repeat filter (default: repfiltered/fastq.gz.stats)
#  OUTDIR= Directory to output files (default: BWAPATH)
#  OUTSTUB= Suffix for all output files (default: counts.txt)
#  CORES= The number of CPUs to use for parallelization (defaults to 1)
#   NOTE: If submitting your job through SLURM, make sure to set -c option as well
#
# UPDATE 3/4/16 - Collapsing counts from re-sequenced library in this script
# NOTE: This means the script must now be run AFTER build_sample_table.R
#   The script does NOT need to be given the batch list anymore,
#   that information is now pulled from [library/sample_master]_table.txt files
#

# setwd("~/Projects/DGRP_Baseline_RNAseq_Align")

# For melt functions:
library(reshape2)
# For parallelization:
library(doMC)

# library(R.utils)

options(stringsAsFactors=F)

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript repeat_counts.R [OPTIONS]"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("GROUPS=~/Resources/RepBase/dmel.repbase-20.01.annot.txt")
# my.args <- c("GROUPS=~/Resources/RepBase/dmel.repbase-20.01.annot.txt")
# my.args <- c("GROUPS=~/Resources/RepBase/dmel.repbase-20.01.annot.txt", "SAMPLES=DNA_sample_master_table.txt", "LIBRARIES=DNA_library_master_table.txt", "PRE=microbefiltered/DGRP_line_DNA/fastq.gz.stats", "POST=repfiltered/DGRP_line_DNA/fastq.gz.stats", "OUTSTUB=DNA_counts.txt")
if(length(my.args) > 0) {
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

setDefault("SAMPLES","sample_master_table.txt")

setDefault("LIBRARIES","library_master_table.txt")

setDefault("GROUPS","/home/ljeveret/Resources/RepBase/dmel.repbase-20.01.annot.txt")

setDefault("BWAPATH","bwa_repeat/")
# Add trailing / if missing
if(!grepl("[/]$", my.args$BWAPATH)) my.args$BWAPATH <- paste0(my.args$BWAPATH, "/")
# Make sure this subdir exists
if(!file.exists(my.args$BWAPATH)) {
  stop("BWAPATH = ", my.args$BWAPATH, " does not exist! Have you run BWA against RepBase yet?")
}

setDefault("COUNTS","target_counts.txt")

# OUTDIR default: Extract path to EXPR input file, use it as the output dir, append Perm if doing permutation analysis
setDefault("OUTDIR", my.args$BWAPATH)
if(!grepl("[/]$", my.args$OUTDIR)) my.args$OUTDIR <- paste0(my.args$OUTDIR, "/")
# Make sure this subdir exists - create if it doesn't
if(!file.exists(my.args$OUTDIR)) {
  dir.create(path=my.args$OUTDIR, recursive=T, showWarnings = F)
}

setDefault("PRE", "microbefiltered/fastq.gz.stats")

setDefault("POST", "repfiltered/fastq.gz.stats")

setDefault("OUTSTUB", "counts.txt")

# Set number of cores for multi-tasking
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
  my.args$CORES <- as.integer(1)
  cat("Setting CORES = 1 by default.\n")
}
# If CORES = ALL, use detectCores()-1
if(is.character(my.args$CORE) & (my.args$CORE == "ALL")) {
  my.args$CORES <- as.integer(detectCores()-1)
  cat("Setting CORES =", my.args$CORES, "[ALL]\n")
}
# Otherwise, attempt to convert to integer
if(!is.integer(my.args$CORE)) {
  my.args$CORE <- as.integer(my.args$CORE)
  cat("Setting CORES =", my.args$CORE, "\n")
}

# Set up parallel backend with desired number of CPUs
cat("Using", my.args$CORES, "CPUs for parallel tasks.\n")
registerDoMC(my.args$CORES)
cat("\n\n")


# --- LOAD LIBRARY AND SAMPLE TABLES --- #

# Load the library table - has batch data, quit if it can't be found
if(!file.exists(my.args$LIBRARIES)) {
  stop("Missing ", my.args$LIBRARIES, " - run build_sample_table.R first!\n")
}
cat("Loading library info:",my.args$LIBRARIES,"\n")
library.table <- read.table(my.args$LIBRARIES, header=T, sep="\t", row.names=1)

# Load the sample table - needed to aggregate counts for re-sequenced libraries
if(!file.exists(my.args$SAMPLES)) {
  stop("Missing ", my.args$SAMPLES, " - run build_sample_table.R first!\n")
}
cat("Loading sample annotation info:",my.args$SAMPLES,"\n")
sample.table <- read.table(my.args$SAMPLES, header=T, sep="\t", row.names=1)
sample.library.map <- strsplit(as.character(sample.table$LIBRARIES), split=",", fixed=T)
names(sample.library.map) <- row.names(sample.table)
cat("\n")

# Check for other files (load them later)
if(!file.exists(my.args$PRE)) {
  stop("Missing ", my.args$PRE, " - run file_stats.sh on microbefiltered fastq files first!")
}
if(!file.exists(my.args$POST)) {
  stop("Missing ", my.args$POST, " - run file_stats.sh on repfiltered fastq files first!")
}


# --- META-INFO FOR REPBASE SEQUENCES --- #

# Load class info for summarizing across target sequences
cat("Loading RepBase sequence class info from:\n")
cat(my.args$GROUPS, "\n")
repeat.info <- read.table(my.args$GROUPS, sep="\t", row.names=1)
cat("Loaded information for",nrow(repeat.info),"primary sequence IDs.\n")

# -- Build mapping data structures for summarizing read groups -- #
# Sequence ID -> Class map
seq.class.map <- repeat.info[,1]
names(seq.class.map) <- row.names(repeat.info)

# Function for applying any of these maps to ID sets
group.mapper <- function(ids, map, delim=","){
  id.list <- strsplit(ids, delim, fixed=T)
  id.list <- lapply(id.list, function(x){sort(unique(map[x]))})
  return(unlist(lapply(id.list, paste, collapse=",")))
}

seq.class.mapper <- function(ids, delim=","){group.mapper(ids, map=seq.class.map, delim=delim)}


# --- LOAD LIBRARY REPBASE COUNT DATA --- #

# Loop over libraries, load the repeat count table for each one, stitch them together into a single table
# Note every sample has the same rows, so need to be careful here
# TO DO: This code can maybe be sped up by loading into a list structure, then building the union of unique IDs?
libs <- row.names(library.table)
cat("Loading target-level count data for",length(libs),"libraries.\n")
# Load all count data into list structure first (because not all target IDs will match up, merge function is too slow here)
count.list <- list()
for(lib in libs) {
  lib.col <- read.table(paste0(my.args$BWAPATH,library.table[lib,"BATCH"],"/",lib,"/",my.args$COUNTS), sep="\t")
  count.list[[lib]] <- lib.col[,2]
  names(count.list[[lib]]) <- lib.col[,1]
}
# Extract the list of all unique target ID sets
target.IDs <- sort(unique(unlist(lapply(count.list, names))))
# Now construct a final data frame
count.table <- data.frame(row.names=target.IDs)
for(lib in libs) {
  lib.targets <- names(count.list[[lib]])
  count.table[lib.targets,lib] <- count.list[[lib]]
}
# Make sure all IDs are in repeat.info table
row.ID.split <- strsplit(row.names(count.table), ",", fixed=T)
unique.IDs <- unique(c(row.ID.split, recursive=T))
stopifnot(all(unique.IDs %in% row.names(repeat.info)))
# Replace NA values with 0
for(j in 1:ncol(count.table)) {
  count.table[is.na(count.table[,j]),j] <- 0
}
for(j in 1:ncol(count.table)) {
  stopifnot(sum(is.na(count.table[,j]))==0)
}

# Check library sums against total reads computed from fastq files
# Load the fastq.gz.wcl tables for both ribofiltered and microbefiltered
microbe.fastq.wcl <- read.table(my.args$PRE, sep="\t", row.names=1)
# Parse out lib name for RNA-seq samples
row.names(microbe.fastq.wcl) <- sub("^.*/","",sub("_filtered.fastq.gz","",row.names(microbe.fastq.wcl)))
# Parse out lib name for DNA-seq samples
row.names(microbe.fastq.wcl) <- sub("^.*/","",sub("_DNA_orig.fastq.gz","",row.names(microbe.fastq.wcl)))
if(!all(libs %in% row.names(microbe.fastq.wcl))) {
  missing.libs <- setdiff(libs, row.names(microbe.fastq.wcl))
  cat(length(missing.libs), "of", length(libs), "library names don't match up to", nrow(microbe.fastq.wcl), "microbe fastq files.\n")
  cat("Library names look like:\n")
  print(head(missing.libs))
  unmatched.fastq <- setdiff(row.names(microbe.fastq.wcl), libs)
  cat("Unmatched Microbe fastq file names parsed to:\n")
  print(head(unmatched.fastq))
  quit(save="no")
}

repeat.fastq.wcl <- read.table(my.args$POST, sep="\t", row.names=1)
row.names(repeat.fastq.wcl) <- sub("^.*/","",sub("_filtered.fastq.gz","",row.names(repeat.fastq.wcl)))
row.names(repeat.fastq.wcl) <- sub("^.*/","",sub("_DNA_orig.fastq.gz","",row.names(repeat.fastq.wcl)))
if(!all(libs %in% row.names(repeat.fastq.wcl))) {
  missing.libs <- setdiff(libs, row.names(repeat.fastq.wcl))
  cat(length(missing.libs), "of", length(libs), "library names don't match up to", nrow(repeat.fastq.wcl), "repeat fastq files.\n")
  cat("Library names look like:\n")
  print(head(missing.libs))
  unmatched.fastq <- setdiff(row.names(repeat.fastq.wcl), libs)
  cat("Unmatched repeat fastq file names parsed to:\n")
  print(head(unmatched.fastq))
  quit(save="no")
}

# compute number of microbe-matching reads
repeat.total.reads <- as.integer((microbe.fastq.wcl[libs,1] - repeat.fastq.wcl[libs,1]) / 4)
names(repeat.total.reads) <- libs

# Compare to column sums in repeat read count table
cat("\n")
count.table.totals <- as.integer(round(apply(count.table, 2, sum)))
names(count.table.totals) <- libs
if(all(count.table.totals == repeat.total.reads)) {
  cat("Repeat counts match up to difference between PRE and POST filters exactly - these are unfiltered counts.")
} else if (all(count.table.totals <= repeat.total.reads)) {
  # Compute % filtered for each library, report the distribution of those values
  perc.filtered <- (repeat.total.reads - count.table.totals) / repeat.total.reads
  cat("Some repeat alignments had reads filtered (could be due to paired end reads with same name).\nDistribution of filter rates:\n")
  cat(min(perc.filtered), boxplot.stats(perc.filtered)$stats, max(perc.filtered), "\n")
} else {
  # Error condition - some microbe counts add up to more than max possible
  bad.libs <- names(count.table.totals)[which(count.table.totals > repeat.total.reads)]
  stop(length(bad.libs), " libraries have more total reads than indicated by PRE and POST fastq counts:\n", paste(bad.libs, collapse="\n"))
}

cat("\nAll count data loaded and validated.\n\n")


# --- Total counts across all librarie for each repbase target sequence --- #

cat("Summarizing total read counts across all",length(libs),"libraries.\n")

# Create mapping table for sequence *sets*
target.totals <- data.frame(CLASS=seq.class.mapper(row.names(count.table)), row.names=row.names(count.table))
stopifnot(sum(is.na(target.totals$CLASS))==0)

# Sum across all libs (no weighting by total sequences, just want to know what the top grossing sequences are)
target.totals[,"TOTAL"] <- apply(count.table, 1, sum)

# What are the top sequence sets by total count?
cat("Top individual target sequence sets:\n")
print(head(target.totals[order(target.totals$TOTAL, decreasing=T),]))
cat("\n")
# Based on all batches:
# COPIA_DM_I still dominates, ROO_I is 2nd with similar order

# Are any RepBase sequences completely missing?
all.targets <- unique(c(strsplit(row.names(target.totals), split=",", fixed=T), recursive=T))
missing.targets <- setdiff(row.names(repeat.info), all.targets)
if(length(missing.targets) > 0) {
  cat("No alignments to", length(missing.targets), "sequences in RepBase:\n")
  print(missing.targets)
  cat("\n")
}

# Did any RepBase targets fail to align uniquely to any reads?
unique.targets <- grep(",", row.names(target.totals), fixed=T, invert=T, value=T)
ambig.targets <- setdiff(row.names(repeat.info),unique.targets)
if(length(ambig.targets) > 0) {
  cat("No unique alignments to", length(ambig.targets), "sequences in RepBase:\n")
  print(ambig.targets)
  cat("\n")
}

# What percent of total reads are ambiguously aligned to multiple targets?
ambig.perc <- sum(target.totals[grepl(",",row.names(target.totals)),"TOTAL"]) * 100 / sum(target.totals$TOTAL)
cat(ambig.perc, "% of reads align ambiguously to multiple target sequences.\n\n", sep="")

# Write out to disk
total.target.file <- paste0(my.args$OUTDIR, "total_target_", my.args$OUTSTUB)
target.total.output <- cbind(TARGET=row.names(target.totals), target.totals)
target.total.output <- target.total.output[order(target.total.output$TOTAL, decreasing=T),]
cat("Writing target total counts to:", total.target.file, "\n")
if(file.exists(total.target.file)) {
  cat("WARNING:", total.target.file, "exists already - Overwriting!\n")
}
write.table(target.total.output, total.target.file, sep="\t", row.names=F, quote=F)

# Summarize by class (use melt and dcast)
cat("Summarizing total counts by repeat class.\n")
class.melt <- melt(target.totals[,c("CLASS","TOTAL")], id.vars="CLASS")
class.total.counts <- dcast(class.melt, CLASS ~ variable, fun.aggregate=sum)

# Show the top classes by total count
cat("Top RepBase hit classes:\n")
print(head(class.total.counts[order(class.total.counts$TOTAL, decreasing=T),]))
cat("\n")
# Based on all batches,
# Top hits is Gypsy, then Copia, BEL, Jockey (drops off ~5x after that)

# What % are ambiguous at class level?
ambig.perc <- sum(class.total.counts[grepl(",",class.total.counts$CLASS),"TOTAL"]) * 100 / sum(class.total.counts$TOTAL)
cat(round(ambig.perc, digits=1), "% of reads align ambiguously to multiple target sequences.\n\n", sep="")


# Output the class total counts
total.class.file <- paste0(my.args$OUTDIR, "total_class_", my.args$OUTSTUB)
cat("Writing class total counts to:", total.class.file, "\n")
if(file.exists(total.class.file)) {
  cat("WARNING:", total.class.file, "exists already - Overwriting!\n")
}
write.table(class.total.counts[order(class.total.counts$TOTAL, decreasing=T),], total.class.file, sep="\t", row.names=F, quote=F)


# --- Individual Library Analysis --- #

cat("\nAnalyzing counts for individual libraries.\n\n")

# Test for targets w/ < 100 UNIQUE reads in ALL libraries
target.max <- apply(count.table, 1, max)
low.targets <- names(target.max)[target.max < 100]
# Filter to unique assembly hits
low.targets <- grep(",", low.targets, fixed=T, invert=T, value=T)
cat(length(low.targets), "targets have < 100 unique reads in ALL libraries.\n")
print(target.totals[low.targets,])
cat("\n")

# Write this table out to file:
count.table.sums <- apply(count.table, 1, sum)
count.table.output <- cbind(TARGET=row.names(count.table), count.table)
count.table.output <- count.table.output[order(count.table.sums, decreasing=T),]
target.count.file <- paste0(my.args$BWAPATH, "library_target_", my.args$OUTSTUB)
cat("Writing library target-level counts to:", target.count.file, "\n")
if(file.exists(target.count.file)) {
  cat("WARNING:", target.count.file, "already exists - Overwriting!\n")
}
write.table(count.table.output, target.count.file, sep="\t", row.names=F, quote=F)

# Aggregate counts for re-sequenced libraries
cat("Aggregating target counts for each unique sample.\n")
sample.target.count.table <- data.frame(row.names=row.names(count.table))
library.corr.table <- data.frame(SAMPLE=character(0), LIB1=character(0), LIB2=character(0), SCC=numeric(0))
for(smp in names(sample.library.map)) {
  if(length(sample.library.map[[smp]])==1) {
    # When a sample was sequenced only once, just copy the column over to this table
    sample.target.count.table <- cbind(sample.target.count.table, count.table[,sample.library.map[[smp]]])
  } else {
    # When a sample has multiple libraries, check all pairwise correlations and then sum the counts across libraries
    smp.lib.counts <- count.table[,sample.library.map[[smp]]]
    # If more than 1 column, check all parwise correlations (spearman)
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
    sample.target.count.table <- cbind(sample.target.count.table, apply(smp.lib.counts, 1, sum))
    rm(smp.lib.counts)
  }  
}
colnames(sample.target.count.table) <- names(sample.library.map)

# Check the library-wise correlations - as long as these are all very high, there is no need to look into library-level batch effects
# NOTE Pearson Correlation is quite high, but Spearman Correlation levels are low
# This is probably because most target sequences are very low, essentially noise, so I am ignoring the SCC levels here
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

# What % of reads are ambiguously lined per sample on avg
sample.totals <- apply(sample.target.count.table, 2, sum)
sample.ambig <- apply(sample.target.count.table[grepl(",",row.names(sample.target.count.table)),], 2, sum)
sample.ambig.perc <- sample.ambig * 100 / sample.totals
cat(round(mean(sample.ambig.perc)), "% of reads per sample are ambiguous at target level on average.\n\n", sep="")

# Output the sample-level counts
sample.table.sums <- apply(sample.target.count.table, 1, sum)
sample.target.count.table <- cbind(TARGET=row.names(sample.target.count.table), as.data.frame(sample.target.count.table))
sample.target.count.table <- sample.target.count.table[order(sample.table.sums, decreasing=T),]
target.sample.file <- paste0(my.args$BWAPATH, "sample_target_", my.args$OUTSTUB)
cat("Writing sample target-level counts to:", target.sample.file, "\n")
if(file.exists(target.sample.file)) {
  cat("WARNING:", target.sample.file, "already exists - Overwriting!\n")
}
write.table(sample.target.count.table, target.sample.file, sep="\t", row.names=F, quote=F)


# --- TARGET CLASS-LEVEL ANALYSIS --- #

# Now keep the counts separated by library, but condense by class
stopifnot(all(row.names(count.table) %in% row.names(target.totals)))
count.table <- cbind(CLASS=target.totals[row.names(count.table),"CLASS"], count.table)

# Collapse by class
cat("Collapsing library counts by repeat class.\n")
class.melt <- melt(count.table, id.vars="CLASS")
class.count.table <- dcast(class.melt, CLASS ~ variable, fun.aggregate=sum)
row.names(class.count.table) <- class.count.table[,"CLASS"]
class.count.table <- class.count.table[,-1]

# For each library, compute the most abundant class set
library.top.class <- apply(class.count.table, 2, function(x){row.names(class.count.table)[which.max(x)]})

# Add the count of libraries in which each sequence is the top hit
class.total.counts[,"TOP.LIB.COUNT"] <- unlist(lapply(class.total.counts$CLASS, function(x){sum(library.top.class==x)}))
top.lib.n <- sum(class.total.counts$TOP.LIB.COUNT > 0)
cat(top.lib.n,"class sets are the top hit in at least one library:\n")
print(class.total.counts[order(class.total.counts$TOP.LIB.COUNT, decreasing=T)[1:top.lib.n],])
cat("\n")

# Test for classes w/ < 1,000 UNIQUE reads in ALL libraries
class.max <- apply(class.count.table, 1, max)
low.class <- names(class.max)[class.max < 1000]
# Filter to unique class hits
low.class <- grep(",", low.class, fixed=T, invert=T, value=T)
cat(length(low.class), "repeat classes have < 1,000 unique reads in ALL libraries.\n")
print(class.total.counts[class.total.counts$CLASS %in% low.class,])
cat("\n")

# OUTPUT THE RAW COUNTS, SIMILAR FORMAT TO GENE READ COUNTS
# NOTE: It might make more sense to look at reads as proportion of total reads fed into this alignment (ribofiltered total)

# Reorder by sum
stopifnot(row.names(class.count.table) == class.total.counts$CLASS)
class.count.output <- cbind(CLASS=row.names(class.count.table), class.count.table)
class.count.output <- class.count.output[order(class.total.counts$TOTAL, decreasing=T),]

# Write out
library.class.file <- paste0(my.args$OUTDIR, "library_class_", my.args$OUTSTUB)
cat("Writing library class-level counts to:", library.class.file, "\n")
if(file.exists(library.class.file)) {
  cat("WARNING:", library.class.file, "exists already - Overwriting!\n")
}
write.table(class.count.output, library.class.file, sep="\t", row.names=F, quote=F)

# Aggregate counts for re-sequenced libraries
cat("Aggregating class counts for each unique sample.\n")
sample.class.count.table <- data.frame(row.names=row.names(class.count.table))
library.corr.table <- data.frame(SAMPLE=character(0), LIB1=character(0), LIB2=character(0), SCC=numeric(0))
for(smp in names(sample.library.map)) {
  if(length(sample.library.map[[smp]])==1) {
    # When a sample was sequenced only once, just copy the column over to this table
    sample.class.count.table <- cbind(sample.class.count.table, class.count.table[,sample.library.map[[smp]]])
  } else {
    # When a sample has multiple libraries, check all pairwise correlations and then sum the counts across libraries
    smp.lib.counts <- class.count.table[,sample.library.map[[smp]]]
    # If more than 1 column, check all parwise correlations (spearman)
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
    sample.class.count.table <- cbind(sample.class.count.table, apply(smp.lib.counts, 1, sum))
    rm(smp.lib.counts)
  }  
}
colnames(sample.class.count.table) <- names(sample.library.map)

# Check the library-wise correlations - as long as these are all very high, there is no need to look into library-level batch effects
# NOTE Pearson Correlation is quite high, but Spearman Correlation levels are low
# This is probably because most target sequences are very low, essentially noise, so I am ignoring the SCC levels here
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

# What % of reads are ambiguously aligned per sample on avg
sample.totals <- apply(sample.class.count.table, 2, sum)
sample.ambig <- apply(sample.class.count.table[grepl(",",row.names(sample.class.count.table)),], 2, sum)
sample.ambig.perc <- sample.ambig * 100 / sample.totals
cat(round(mean(sample.ambig.perc), digits=1), "% of reads per sample are ambiguous at class level on average.\n\n", sep="")

# Output the sample-level counts
sample.table.sums <- apply(sample.class.count.table, 1, sum)
sample.class.count.table <- cbind(TARGET=row.names(sample.class.count.table), as.data.frame(sample.class.count.table))
sample.class.count.table <- sample.class.count.table[order(sample.table.sums, decreasing=T),]
class.sample.file <- paste0(my.args$BWAPATH, "sample_class_", my.args$OUTSTUB)
cat("Writing sample class-level counts to:", class.sample.file, "\n")
if(file.exists(class.sample.file)) {
  cat("WARNING:", class.sample.file, "already exists - Overwriting!\n")
}
write.table(sample.class.count.table, class.sample.file, sep="\t", row.names=F, quote=F)

