#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 12/16/15
# 
# GOAL: Collect the read counts from BWA alignment against DGRP_Microbiome
# Summarize the counts across all samples and within samples (percentage) and within taxa
#
# Usage:
# Rscript microbe_counts.R [OPTIONS]
#  SAMPLES= The table with sample information (default: sample_master_table.txt)
#  LIBRARIES= The table with library batch information (default: library_master_table.txt)
#  GROUPS= The table of multiple grouping levels (assembly, species, etc.) for all target sequences in microbe DB
#    Default: /home/ljeveret/Resources/NCBI/DGRP_Microbiome/DGRP_Microbiome.groups.txt
#  BWAPATH= The path to top-level BWA results where the count files live (default: bwa_microbe/)
#  COUNTS= The name of count files in each BWA output directory (default: target_counts.txt)
#  PRE= The file containing fastq read counts BEFORE microbe filter (default: ribofiltered/fastq.gz.stats)
#  POST= The file containing fastq read counts AFTER microbe filter (default: microbefiltered/fastq.gz.stats)
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
usageStr="USAGE:\nRscript microbe_counts.R [OPTIONS]"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("GROUPS=~/Resources/NCBI/DGRP_Microbiome/DGRP_Microbiome.groups.txt")
# my.args <- c("GROUPS=~/Resources/NCBI/DGRP_Microbiome/DGRP_Microbiome.groups.txt", "COUNTS=target_filtered_counts.txt", "OUTSTUB=filtered_counts.txt", "CORES=4")
# my.args <- c("GROUPS=~/Resources/NCBI/DGRP_Microbiome/DGRP_Microbiome.groups.txt", "SAMPLES=DNA_sample_master_table.txt", "LIBRARIES=DNA_library_master_table.txt", "COUNTS=target_filtered_counts.txt", "PRE=dna/fastq.gz.stats", "POST=microbefiltered/DGRP_line_DNA/fastq.gz.stats", "OUTSTUB=DNA_counts.txt")
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

setDefault("GROUPS","/home/ljeveret/Resources/NCBI/DGRP_Microbiome/DGRP_Microbiome.groups.txt")

setDefault("BWAPATH","bwa_microbe/")
# Add trailing / if missing
if(!grepl("[/]$", my.args$BWAPATH)) my.args$BWAPATH <- paste0(my.args$BWAPATH, "/")
# Make sure this subdir exists
if(!file.exists(my.args$BWAPATH)) {
  stop("BWAPATH = ", my.args$BWAPATH, " does not exist! Have you run BWA against microbe DB yet?")
}

setDefault("COUNTS","target_counts.txt")

# OUTDIR default: Extract path to EXPR input file, use it as the output dir, append Perm if doing permutation analysis
setDefault("OUTDIR", my.args$BWAPATH)
if(!grepl("[/]$", my.args$OUTDIR)) my.args$OUTDIR <- paste0(my.args$OUTDIR, "/")
# Make sure this subdir exists - create if it doesn't
if(!file.exists(my.args$OUTDIR)) {
  dir.create(path=my.args$OUTDIR, recursive=T, showWarnings = F)
}

setDefault("PRE", "ribofiltered/fastq.gz.stats")

setDefault("POST", "microbefiltered/fastq.gz.stats")

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
  stop("Missing ", my.args$PRE, " - run file_stats.sh on ribofiltered fastq files first!")
}
if(!file.exists(my.args$POST)) {
  stop("Missing ", my.args$POST, " - run file_stats.sh on microbefiltered fastq files first!")
}


# --- META-INFO FOR MICROBE DB --- #

# Load group info for summarizing across assemblies, species, etc.
cat("Loading DGRP_Microbiome sequence and grouping info from:\n")
cat(my.args$GROUPS, "\n")
group.info <- read.table(my.args$GROUPS, header=T, sep="\t", row.names=1)
cat("Loaded information for",nrow(group.info),"target sequence IDs.\n")

# -- Build mapping data structures for summarizing read groups -- #
cat("Constructing meta-group mappings...\n")

# Refseq -> Assembly map
refseq.assembly.map <- group.info$ASSEMBLY
names(refseq.assembly.map) <- row.names(group.info)

# Assembly -> Organism map
assembly.org.pairs <- paste(group.info$ASSEMBLY, group.info$ORGANISM, sep=":")
assembly.org.unique <- which(!duplicated(assembly.org.pairs))
stopifnot(length(assembly.org.unique)==length(unique(group.info$ASSEMBLY)))
assembly.org.map <- group.info[assembly.org.unique,"ORGANISM"]
names(assembly.org.map) <- group.info[assembly.org.unique,"ASSEMBLY"]
rm(assembly.org.pairs, assembly.org.unique)
# Make sure the current group separator does not conflict with any organism names
stopifnot(sum(grepl(",",assembly.org.map))==0)

# Assembly -> TaxID map
assembly.taxid.pairs <- paste(group.info$ASSEMBLY, group.info$TAXID, sep=":")
assembly.taxid.unique <- which(!duplicated(assembly.taxid.pairs))
stopifnot(length(assembly.taxid.unique)==length(unique(group.info$ASSEMBLY)))
assembly.taxid.map <- group.info[assembly.taxid.unique,"TAXID"]
names(assembly.taxid.map) <- group.info[assembly.taxid.unique,"ASSEMBLY"]
rm(assembly.taxid.pairs, assembly.taxid.unique)

# Organism -> Species map
org.species.pairs <- paste(group.info$ORGANISM, group.info$SPECIES, sep=":")
org.species.unique <- which(!duplicated(org.species.pairs))
stopifnot(length(org.species.unique)==length(unique(group.info$ORGANISM)))
org.species.map <- group.info[org.species.unique,"SPECIES"]
names(org.species.map) <- group.info[org.species.unique,"ORGANISM"]
rm(org.species.pairs, org.species.unique)
# Make sure the current group separator does not conflict with any species names
stopifnot(sum(grepl(",",org.species.map))==0)

# Species -> Genus map
species.genus.pairs <- paste(group.info$SPECIES, group.info$GENUS, sep=":")
species.genus.unique <- which(!duplicated(species.genus.pairs))
stopifnot(length(species.genus.unique)==length(unique(group.info$SPECIES)))
species.genus.map <- group.info[species.genus.unique,"GENUS"]
names(species.genus.map) <- group.info[species.genus.unique,"SPECIES"]
rm(species.genus.pairs, species.genus.unique)
stopifnot(sum(grepl(",",species.genus.map))==0)

# Function for applying any of these maps to ID sets
group.mapper <- function(ids, map, delim=","){
  id.list <- strsplit(ids, delim, fixed=T)
  id.list <- lapply(id.list, function(x){sort(unique(map[x]))})
  return(unlist(lapply(id.list, paste, collapse=",")))
}
cat("Group maps complete.\n\n")


# --- LOAD LIBRARY COUNT DATA --- #

# Loop over samples, load the microbial count table for each one, stitch them together into a single table
# Note every sample has the same rows, so need to be careful here
# TO DO: This code can maybe be sped up by loading into a list structure, then building the union of unique IDs?
libs <- row.names(library.table)
cat("Computing assembly-level count data for",length(libs),"libraries.\n")
# Load all count data into list structure first (because not all target IDs will match up, merge function is too slow here)
# UPDATE 10/19/16 - Counts are immediately collapsed to assembly level to conserve memory and speed up the table merge step
count.list <- foreach(lib=libs) %dopar% {
  count.file <- paste0(my.args$BWAPATH,library.table[lib,"BATCH"],"/",lib,"/",my.args$COUNTS)
  if(!file.exists(count.file)) {
    stop("Missing required count file: ", count.file, " - have you run the BWA target count yet?")
  }
  lib.col <- read.table(count.file, sep="\t", row.names=1)
  
  # Remove gi portion from row names
  # Now need to do this recursively for group names
  # First, split into separate IDs on commas
  row.ID.split <- strsplit(row.names(lib.col), ",", fixed=T)
  # Strip out the leading gi portion and all "|" dividers
  row.ID.split <- lapply(row.ID.split, function(x){
    sub("^gi[|][0-9]*[|]ref[|]","",sub("[|]$","",x))
  })
  # Make sure all IDs are in group.info table
  unique.IDs <- c(row.ID.split, recursive=T)
  stopifnot(all(unique.IDs %in% row.names(group.info)))
  # Paste back together, keep comma as delimiter, use as row names
  refseq.IDs <- unlist(lapply(row.ID.split, paste, collapse=","))
  stopifnot(length(refseq.IDs) == nrow(lib.col))
  stopifnot(all(refseq.IDs != ""))
  stopifnot(all(!is.na(refseq.IDs)))
  row.names(lib.col) <- refseq.IDs
  
  # Now map to assembly and melt the counts
  lib.col$ASSEMBLY <- group.mapper(row.names(lib.col), refseq.assembly.map)
  lib.melt <- melt(lib.col, id.vars="ASSEMBLY")
  lib.col <- dcast(lib.melt, ASSEMBLY ~ variable, fun.aggregate=sum)
  lib.counts <- lib.col[,2]
  names(lib.counts) <- lib.col[,1]
  return(lib.counts)
}
names(count.list) <- libs
# rm(lib.melt, lib.col)

# Extract the list of all unique target IDs
target.IDs <- sort(unique(unlist(lapply(count.list, names))))
# Now construct a final data frame
count.table <- data.frame(row.names=target.IDs)
for(lib in libs) {
  lib.targets <- names(count.list[[lib]])
  count.table[lib.targets,lib] <- count.list[[lib]]
}
# Replace NA values with 0
for(j in 1:ncol(count.table)) {
  count.table[is.na(count.table[,j]),j] <- 0
}
for(j in 1:ncol(count.table)) {
  stopifnot(sum(is.na(count.table[,j]))==0)
}

# Memory clean-up
# rm(count.list, refseq.IDs, row.ID.split, target.IDs)
rm(count.list, target.IDs)

# TO DO: It would make more sense to have these values compiled in library_master_table.txt ahead of time...

# Check library sums against total reads computed from fastq files
# Load the fastq.gz.stats tables for both ribofiltered and microbefiltered
ribo.fastq.wcl <- read.table(my.args$PRE, sep="\t", row.names=1)
# Parse out lib name for RNA-seq samples
row.names(ribo.fastq.wcl) <- sub("^.*/","",sub("_filtered.fastq.gz","",row.names(ribo.fastq.wcl)))
# Parse out lib name for DNA-seq samples
row.names(ribo.fastq.wcl) <- sub("^.*/","",sub("_DNA_orig.fastq.gz","",row.names(ribo.fastq.wcl)))
if(!all(libs %in% row.names(ribo.fastq.wcl))) {
  missing.libs <- setdiff(libs, row.names(ribo.fastq.wcl))
  cat(length(missing.libs), "of", length(libs), "library names don't match up to", nrow(ribo.fastq.wcl), "rRNA fastq files.\n")
  cat("Library names look like:\n")
  print(head(missing.libs))
  unmatched.fastq <- setdiff(row.names(ribo.fastq.wcl), libs)
  cat("Unmatched rRNA fastq file names parsed to:\n")
  print(head(unmatched.fastq))
  quit(save="no")
}

microbe.fastq.wcl <- read.table(my.args$POST, sep="\t", row.names=1)
row.names(microbe.fastq.wcl) <- sub("^.*/","",sub("_filtered.fastq.gz","",row.names(microbe.fastq.wcl)))
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

# compute number of microbe-matching reads
microbe.total.reads <- as.integer((ribo.fastq.wcl[libs,1] - microbe.fastq.wcl[libs,1]) / 4)
names(microbe.total.reads) <- libs

# Compare to column sums in microbe read count table
cat("\n")
count.table.totals <- as.integer(round(apply(count.table, 2, sum)))
names(count.table.totals) <- libs
if(all(count.table.totals == microbe.total.reads)) {
  cat("Microbe counts match up to difference between PRE and POST filters exactly - these are unfiltered counts.")
} else if (all(count.table.totals <= microbe.total.reads)) {
  # Compute % filtered for each library, report the distribution of those values
  perc.filtered <- (microbe.total.reads - count.table.totals) / microbe.total.reads
  cat("Some microbe alignments had reads filtered.\nDistribution of filter rates:\n")
  cat(min(perc.filtered), boxplot.stats(perc.filtered)$stats, max(perc.filtered), "\n")
} else {
  # Error condition - some microbe counts add up to more than max possible
  bad.libs <- names(count.table.totals)[which(count.table.totals > microbe.total.reads)]
  stop(length(bad.libs), " libraries have more total reads than indicated by PRE and POST fastq counts:\n", paste(bad.libs, collapse="\n"))
}

cat("\nAll count data loaded and validated.\n\n")

# Write the assembly by library table out to disk
lib.assembly.file <- paste0(my.args$OUTDIR, "library_assembly_", my.args$OUTSTUB)
cat("Writing library assembly-level count table to:", lib.assembly.file, "\n")
if(file.exists(lib.assembly.file)) {
  cat("WARNING:", lib.assembly.file, "already exists, Overwriting!")
}
write.table(cbind(ASSEMBLY=row.names(count.table), count.table), lib.assembly.file, sep="\t", row.names=F, quote=F)
cat("\n")


# --- Condense Libraries to Sample-Level Counts --- #

# Aggregate counts for re-sequenced libraries
cat("Aggregating assembly counts for each unique sample.\n")
sample.assembly.count.table <- data.frame(row.names=row.names(count.table))
library.corr.table <- data.frame(SAMPLE=character(0), LIB1=character(0), LIB2=character(0), SCC=numeric(0))
for(smp in names(sample.library.map)) {
  if(length(sample.library.map[[smp]])==1) {
    # When a sample was sequenced only once, just copy the column over to this table
    sample.assembly.count.table <- cbind(sample.assembly.count.table, count.table[,sample.library.map[[smp]]])
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
    sample.assembly.count.table <- cbind(sample.assembly.count.table, apply(smp.lib.counts, 1, sum))
    rm(smp.lib.counts)
  }  
}
colnames(sample.assembly.count.table) <- names(sample.library.map)

# Check the library-wise correlations - as long as these are all very high, there is no need to look into library-level batch effects
# NOTE Pearson Correlation is quite high, but Spearman Correlation levels are low
# This is probably because most target sequences are very low, essentially noise, so I am ignoring the SCC levels here
cat("Checking correlations for re-sequenced libraries:\n\n")
cat("Distribution of Pearson Correlation Coefficients:\n")
cat(min(library.corr.table$PCC), boxplot.stats(library.corr.table$PCC)$stats, max(library.corr.table$PCC), "\n")
cat("Distribution of Spearman Correlation Coefficients:\n")
cat(min(library.corr.table$SCC), boxplot.stats(library.corr.table$SCC)$stats, max(library.corr.table$SCC), "\n")
cat("\n")

# Remove the library count table
rm(count.table)

# Convert sample assembly-level count table to data frame
sample.assembly.count.table <- as.data.frame(sample.assembly.count.table)
samples <- colnames(sample.assembly.count.table)


# --- Total counts across all samples for each microbe assembly --- #

# TO DO: Remove all steps that handle refseq level counts - those are now collapsed to assembly during the individual count file loading step
# TO DO: Should check replicate library correlation at assembly level and then immediately collapse to sample level counts
#       That is, only process/ouput library level count table at assembly level

cat("Summarizing total read counts across all",length(samples),"samples.\n")

# Create mapping table for sequence *sets*
assembly.totals <- data.frame(ASSEMBLY=row.names(sample.assembly.count.table), row.names=row.names(sample.assembly.count.table))
stopifnot(sum(is.na(assembly.totals$ASSEMBLY))==0)
assembly.totals[,"ORGANISM"] <- group.mapper(assembly.totals$ASSEMBLY, assembly.org.map)
stopifnot(sum(is.na(assembly.totals$ORGANISM))==0)
assembly.totals[,"SPECIES"] <- group.mapper(assembly.totals$ORGANISM, org.species.map)
stopifnot(sum(is.na(assembly.totals$SPECIES))==0)
assembly.totals[,"GENUS"] <- group.mapper(assembly.totals$SPECIES, species.genus.map)
stopifnot(sum(is.na(assembly.totals$GENUS))==0)

# Sum across all samples (no weighting by total sequences, just want to know what the top grossing sequences are)
assembly.totals[,"TOTAL"] <- apply(sample.assembly.count.table, 1, sum)

# What are the top sequence sets by total count?
cat("Top assembly sequence sets:\n")
print(head(assembly.totals[order(assembly.totals$TOTAL, decreasing=T),]))
cat("\n")
# Based on all batches:
# Wolbachia dominates, Noravirus is 2nd hit, followed by Aspergillus

# Are any assemblies missing here (no alignments)
all.assemblies <- unique(c(strsplit(assembly.totals$ASSEMBLY, ",", fixed=T), recursive=T))
missing.assemblies <- setdiff(group.info$ASSEMBLY, all.assemblies)
if(length(missing.assemblies) > 0) {
  cat("The following assemblies have no alignments:\n")
  print(group.info[group.info$ASSEMBLY %in% missing.assemblies,])
}

# Do any assemblies fail to get any UNIQUE alignments
# (that is, all their alignments are shared across multiple assemblies)
unique.assemblies <- grep(",", assembly.totals$ASSEMBLY, fixed=T, invert=T, value=T)
ambig.assemblies <- setdiff(all.assemblies,unique.assemblies)
if(length(ambig.assemblies) > 0) {
  cat("The following assemblies only have amiguous alignments:\n")
  print(group.info[group.info$ASSEMBLY %in% ambig.assemblies,])
}

perc.ambig <- sum(assembly.totals[grepl(",",row.names(assembly.totals)),"TOTAL"]) * 100 / sum(assembly.totals$TOTAL)
cat(round(perc.ambig), "% of reads align to multiple genome assemblies.\n\n\n", sep="")

# Summarize by Species
cat("Summarizing total counts by species.\n")
species.melt <- melt(assembly.totals[,c("SPECIES","TOTAL")], id.vars="SPECIES")
species.counts <- dcast(species.melt, SPECIES ~ variable, fun.aggregate=sum)
species.counts <- species.counts[order(species.counts$TOTAL, decreasing=T),]
rm(species.melt)

# Show the top species sets by total count
cat("Top 10 species sets account for",round(sum(species.counts[1:10,"TOTAL"])*100/sum(species.counts$TOTAL)),"% of all reads:\n")
print(species.counts[1:10,])
cat("\n")
# Top species are Wolbachia, Nora Virus, Acetobacter, Aspergillus

perc.ambig <- sum(species.counts[grepl(",",species.counts$SPECIES),"TOTAL"]) * 100 / sum(species.counts$TOTAL)
cat(round(perc.ambig), "% of reads align to multiple species.\n\n\n", sep="")


# Summarize by Genus
cat("Summarizing total counts by genus.\n")
genus.melt <- melt(assembly.totals[,c("GENUS","TOTAL")], id.vars="GENUS")
genus.counts <- dcast(genus.melt, GENUS ~ variable, fun.aggregate=sum)
genus.counts <- genus.counts[order(genus.counts$TOTAL, decreasing=T),]
rm(genus.melt)

# What are top genus groups, and what % do they account for
cat("Top 10 genus sets account for", round(sum(genus.counts[1:10,"TOTAL"])*100/sum(genus.counts$TOTAL)), "% of all reads:\n")
print(genus.counts[1:10,])
cat("\n")
# Based on All Batches:
# Wolbachia tops the list, followed by Acetobacter,
# Then Nora Virus, Aspergillus, Tremella, Podospora, Lactobacillus, Chaemotium, Plasmodium, Pectobacterium

# What % of reads is an ambiguous genus group?
ambig.genus.perc <- sum(genus.counts[grepl(",", genus.counts$GENUS, fixed=T),"TOTAL"])*100/sum(genus.counts$TOTAL)
cat(round(ambig.genus.perc), "% of all reads are ambiguously aligned at genus level.\n")

# Output the genus total counts
total.genus.file <- paste0(my.args$OUTDIR, "total_genus_", my.args$OUTSTUB)
cat("Writing genus total counts to:", total.genus.file, "\n")
if(file.exists(total.genus.file)) {
  cat("WARNING:", total.genus.file, "exists already - Overwriting!\n")
}
write.table(genus.counts, total.genus.file, sep="\t", row.names=F, quote=F)
cat("\n")


# --- Individual Sample Analysis --- #

cat("\nAnalyzing counts for individual samples.\n\n")

# Compute avg % of reads ambiguously aligned per library
sample.totals <- apply(sample.assembly.count.table, 2, sum)
sample.perc.ambig <- apply(sample.assembly.count.table[grepl(",", row.names(sample.assembly.count.table)),], 2, sum)
sample.perc.ambig <- sample.perc.ambig * 100 / sample.totals
cat(round(mean(sample.perc.ambig)), "% of reads per sample align to multiple target sequences on average.\n\n\n", sep="")

# For each sample, compute the most abundant assembly set
sample.top.assembly <- apply(sample.assembly.count.table, 2, function(x){row.names(sample.assembly.count.table)[which.max(x)]})

# Add the count of libraries in which each sequence is the top hit to assembly.counts
assembly.totals[,"TOP.LIB.COUNT"] <- unlist(lapply(assembly.totals$ASSEMBLY, function(x){sum(sample.top.assembly==x)}))
top.lib.n <- sum(assembly.totals$TOP.LIB.COUNT > 0)
cat(top.lib.n,"assembly sets are the top microbe in at least one sample:\n")
print(assembly.totals[order(assembly.totals$TOP.LIB.COUNT, decreasing=T)[1:top.lib.n],])
cat("\n")

# UPDATED 10/20/16: Test for assemblies w/ < 100 UNIQUE reads in ALL libraries
# This threshold is still <1% in ALL samples
assembly.max <- apply(sample.assembly.count.table, 1, max)
low.assembly <- names(assembly.max)[assembly.max < 100]
# Filter to unique assembly hits
low.assembly <- grep(",", low.assembly, fixed=T, invert=T, value=T)
cat(length(low.assembly), "assemblies have < 1,000 unique reads in ALL samples.\n")
print(assembly.totals[assembly.totals$ASSEMBLY %in% low.assembly,c("ASSEMBLY","ORGANISM","TOTAL")])
cat("\n")

# Compute avg % of reads ambiguously aligned per sample
sample.totals <- apply(sample.assembly.count.table, 2, sum)
sample.perc.ambig <- apply(sample.assembly.count.table[grepl(",", row.names(sample.assembly.count.table)),], 2, sum)
sample.perc.ambig <- sample.perc.ambig * 100 / sample.totals
cat(round(mean(sample.perc.ambig)), "% of reads per sample align to multiple genomes on average.\n\n\n", sep="")

# OUTPUT THE RAW COUNTS, SIMILAR FORMAT TO GENE READ COUNTS
# NOTE: It might make more sense to look at reads as proportion of total reads fed into this alignment (ribofiltered total)

# add detail to row names, append as first column
row.IDs <- strsplit(row.names(sample.assembly.count.table), ",", fixed=T)
row.details <- unlist(lapply(row.IDs, function(x){
  paste0(assembly.org.map[x], " [", x, "|", assembly.taxid.map[x],"]", collapse=" // ")
}))
stopifnot(length(row.details)==nrow(sample.assembly.count.table))

assembly.counts.output <- cbind(ASSEMBLY=row.names(sample.assembly.count.table), ORGANISM=row.details, sample.assembly.count.table)

# Reorder by mean %
mean.assembly.reads <- apply(assembly.counts.output[,3:ncol(assembly.counts.output)], 1, mean)
assembly.counts.output <- assembly.counts.output[order(mean.assembly.reads, decreasing=T),]

# Write out
assembly.count.file <- paste0(my.args$BWAPATH, "sample_assembly_", my.args$OUTSTUB)
cat("Writing sample assembly-level counts to:", assembly.count.file, "\n")
if(file.exists(assembly.count.file)) {
  cat("WARNING:", assembly.count.file, "already exists - Overwriting!\n")
}
write.table(assembly.counts.output, assembly.count.file, sep="\t", row.names=F, quote=F)
rm(assembly.counts.output)
cat("\n")


# Collapse by Species
cat("Collapsing sample counts by species.\n")
stopifnot(all(row.names(sample.assembly.count.table) == assembly.totals$ASSEMBLY))
sample.assembly.count.table <- cbind(SPECIES=assembly.totals$SPECIES, sample.assembly.count.table)
species.melt <- melt(sample.assembly.count.table, id.vars="SPECIES")
species.count.table <- dcast(species.melt, SPECIES ~ variable, fun.aggregate=sum)
row.names(species.count.table) <- species.count.table[,"SPECIES"]
species.count.table <- species.count.table[,-1]
rm(species.melt)

# For each sample, compute the most abundant species, count how many times each species is the most abundant
top.species <- table(apply(species.count.table, 2, function(x){row.names(species.count.table)[which.max(x)]}))
cat(length(top.species),"species are the most abundant in at least one sample:\n")
print(top.species[order(top.species, decreasing=T)])
cat("\n")

# UPDATED 10/20/16: Test for assemblies w/ < 100 UNIQUE reads in ALL libraries
# This threshold is still <1% in ALL samples
species.max <- apply(species.count.table, 1, max)
low.species <- names(species.max)[species.max < 100]
low.species <- grep(",", low.species, fixed=T, invert=T, value=T)
cat(length(low.species), "species have < 1,000 unique reads in ALL samples.\n")
print(species.counts[species.counts$SPECIES %in% low.species,c("SPECIES","TOTAL")])
cat("\n")

# Compute avg % of reads ambiguously aligned per sample
sample.totals <- apply(species.count.table, 2, sum)
sample.perc.ambig <- apply(species.count.table[grepl(",", row.names(species.count.table)),], 2, sum)
sample.perc.ambig <- sample.perc.ambig * 100 / sample.totals
cat(round(mean(sample.perc.ambig)), "% of reads per sample align to multiple species on average.\n\n\n", sep="")

# OUTPUT THE RAW COUNTS, SIMILAR FORMAT TO GENE READ COUNTS
# NOTE: It might make more sense to look at reads as proportion of total reads fed into this alignment (ribofiltered total)

# add detail to row names, append as first column
species.counts.output <- cbind(SPECIES=row.names(species.count.table), species.count.table)

# Reorder by mean %
mean.species.reads <- apply(species.counts.output[,2:ncol(species.counts.output)], 1, mean)
species.counts.output <- species.counts.output[order(mean.species.reads, decreasing=T),]

# Write out
species.count.file <- paste0(my.args$BWAPATH, "sample_species_", my.args$OUTSTUB)
cat("Writing sample species-level counts to:", species.count.file, "\n")
if(file.exists(species.count.file)) {
  cat("WARNING:", species.count.file, "already exists - Overwriting!\n")
}
write.table(species.counts.output, species.count.file, sep="\t", row.names=F, quote=F)
rm(species.counts.output, species.count.table)
cat("\n")


# Collapse by genus
cat("Collapsing sample counts by genus.\n")
stopifnot(all(row.names(sample.assembly.count.table) == assembly.totals$ASSEMBLY))
if(colnames(sample.assembly.count.table)[1] == "SPECIES") {
  sample.assembly.count.table <- sample.assembly.count.table[,-1]
}
sample.assembly.count.table <- cbind(GENUS=assembly.totals$GENUS, sample.assembly.count.table)
genus.melt <- melt(sample.assembly.count.table, id.vars="GENUS")
genus.count.table <- dcast(genus.melt, GENUS ~ variable, fun.aggregate=sum)
row.names(genus.count.table) <- genus.count.table[,"GENUS"]
genus.count.table <- genus.count.table[,-1]
rm(genus.melt)

# For each sample, compute the most abundant genus, count how many times each genus is the most abundant
top.genus <- table(apply(genus.count.table, 2, function(x){row.names(genus.count.table)[which.max(x)]}))
cat(length(top.genus),"genus are the most abundant in at least one sample:\n")
print(top.genus[order(top.genus, decreasing=T)])
cat("\n")

# UPDATED 2/13/16: Test for genus w/ < 1,000 UNIQUE reads in ALL libraries
# This threshold is still <1% in ALL samples
genus.max <- apply(genus.count.table, 1, max)
low.genus <- names(genus.max)[genus.max < 1000]
low.genus <- grep(",", low.genus, fixed=T, invert=T, value=T)
cat(length(low.genus), "genera have < 1,000 unique reads in ALL libraries.\n")
print(genus.counts[genus.counts$GENUS %in% low.genus,c("GENUS","TOTAL")])
cat("\n")

# Compute avg % of reads ambiguously aligned per sample
genus.totals <- apply(genus.count.table, 2, sum)
genus.perc.ambig <- apply(genus.count.table[grepl(",", row.names(genus.count.table)),], 2, sum)
genus.perc.ambig <- genus.perc.ambig * 100 / genus.totals
cat(round(mean(genus.perc.ambig)), "% of reads per sample align to multiple genera on average.\n\n\n", sep="")

# OUTPUT THE RAW COUNTS, SIMILAR FORMAT TO GENE READ COUNTS
# NOTE: It might make more sense to look at reads as proportion of total reads fed into this alignment (ribofiltered total)

# add detail to row names, append as first column
genus.counts.output <- cbind(GENUS=row.names(genus.count.table), genus.count.table)

# Reorder by mean %
mean.genus.reads <- apply(genus.counts.output[,2:ncol(genus.counts.output)], 1, mean)
genus.counts.output <- genus.counts.output[order(mean.genus.reads, decreasing=T),]

# Write out
sample.genus.file <- paste0(my.args$BWAPATH, "sample_genus_", my.args$OUTSTUB)
cat("Writing sample genus-level counts to:", sample.genus.file, "\n")
if(file.exists(sample.genus.file)) {
  cat("WARNING:", sample.genus.file, "already exists - Overwriting!\n")
}
write.table(genus.counts.output, sample.genus.file, sep="\t", row.names=F, quote=F)

cat("\n\nScript completed successfully!\n\n")
print(proc.time())
