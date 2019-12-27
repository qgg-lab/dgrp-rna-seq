#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 4/18/16
# Adapted from build_sample_table.R in DGRP_Baseline_Align project
# This version is generalized to handle other project
#
# Build master info table for all samples across all batches
#
# USAGE:
# collect_sample_QC_stats.R BATCH1 [BATCH 2 ...]
#
# TO DO: Need a way to specify a separate table that links sample names to other relevant information (sex, line, treatment)
# These properties cannot be inferred from sample names in generalized way

# This table should contain info for summing read counts from multiple sequencing runs of the same sample
# As well as all batch structure information (including original batch and reseq batch)

options(stringsAsFactors=F)

# Take the list of batches to analyze from command line params
usageStr <- "USAGE:\ncollect_sample_QC_stats.R BATCH1 [BATCH2 ...]"
myArgs <- commandArgs(trailingOnly=T)
if(length(myArgs) == 0) {
  # Can't run without at least one batch!
  stop("Must specify one or more batches to tabulate QC stats. ", usageStr)
} else {
  batches <- myArgs
}

# Remove trailing .txt from these files (just in case)
batches <- sub("[.]txt$", "", batches)

cat("Loading data for",length(batches),"batches:",batches,"\n\n")

# Generate complete sample info table from batch table files
library.table <- NULL
for(batch in batches) {
  batch.file <- paste0(batch,".txt")
  cat("Loading",batch.file,"\n")
  batch.table <- read.table(batch.file, sep="\t")
  # Two different formats allowed for the batch file - ignore the "INDEX" column if missing
  if(ncol(batch.table) == 7) {
    colnames(batch.table) <- c("LIBRARY", "SAMPLE", "FLOWCELL", "LANE", "INDEX", "BARCODE", "SEQDATE")
  } else {
    colnames(batch.table) <- c("LIBRARY", "SAMPLE", "FLOWCELL", "LANE", "BARCODE", "SEQDATE")
  }
  row.names(batch.table) <- batch.table[,"LIBRARY"]
  # Replace "-" in Sample ID with underscore
  batch.table[,"SAMPLE"] <- sub("-", "_", batch.table[,"SAMPLE"])
  # Can't infer these from sample names anymore, need a separate table where these are collected!
  # batch.table[,"LINE"] <- unlist(lapply(strsplit(batch.table[,"SAMPLE"], split="_", fixed=T), function(x){x[1]}))
  # batch.table[,"SEX"] <- as.character(NA)
  # batch.table[grep("_F[1-9]$", batch.table[,"SAMPLE"]),"SEX"] <- "F"
  # batch.table[grep("_M[1-9]$", batch.table[,"SAMPLE"]),"SEX"] <- "M"
  # batch.table[,"REP"] <- as.integer(sub("^[0-9]+_[MF]","",batch.table[,"SAMPLE"]))
  stopifnot(sum(is.na(batch.table)) == 0)
  # stopifnot(all(batch.table[,"REP"] > 0))
  batch.table[,"BATCH"] <- batch
  library.table <- rbind(library.table, batch.table)
}
rm(batch.table)
stopifnot(sum(duplicated(library.table$LIBRARY))==0)
cat("DONE: Loaded information for",nrow(library.table),"individual sequence libraries\n\n")

# Output some basic info: number of flowcells, number of lanes
cat("Complete data set derived from",length(unique(library.table$FLOWCELL)),"flowcells,",length(unique(paste(library.table$FLOWCELL, library.table$LANE))),"individual lanes total.\n\n")

# Drop Index and SeqDate columns...
library.table <- library.table[,setdiff(colnames(library.table),c("INDEX","SEQDATE"))]

# Output this table
library.file <- "library_master_table.txt"
cat("Writing master library file:", library.file, "\n")
write.table(library.table, library.file, row.names=F, sep="\t", quote=F)

# Append some info from sample_QC_stats here: number of input reads at certain alignment steps, number of uniquely aligned reads at the end
lib.QC.table <- read.table("sample_QC_stats.txt", header=T, sep="\t", row.names=1)
stopifnot(all(row.names(library.table) %in% row.names(lib.QC.table)))
lib.QC.table <- lib.QC.table[row.names(library.table),]
stopifnot(all(row.names(lib.QC.table) == row.names(library.table)))
library.table$TOTAL.READS = lib.QC.table[,"TOTAL.READS"]
# Compute the total number of reads fed into microbial alignment (total reads after ribo filtering)
lib.QC.table$TRIMMED.READS = round((1-lib.QC.table$FRAC.READS.TRIMMED)*lib.QC.table$TOTAL.READS)
lib.QC.table$RIBOFILTERED.READS = round((1-lib.QC.table$FRAC.READS.RIBO)*lib.QC.table$TRIMMED.READS)
library.table$MICROBE.LIBSZ = lib.QC.table$RIBOFILTERED.READS
# Compute the total number of reads fed into RepBase alignment (total reads after microbe filtering)
lib.QC.table$MICROBE.FILTERED.READS = round((1-lib.QC.table$FRAC.READS.MICROBE)*lib.QC.table$RIBOFILTERED.READS)
library.table$REPBASE.LIBSZ = lib.QC.table$MICROBE.FILTERED.READS
library.table$LIBSZ = lib.QC.table[,"UNIQUE.ALIGN.READS"]
# TO DO: Could derive some additional technical issue flags here, like high adapter content, high rRNA content, lower alignment rate...

# Now construct sample-level table
samples <- unique(library.table$SAMPLE)
cat("Building sample master info table for",length(samples),"biological samples.\n")

# First construct a list structure, with the list of libraries for each sample
sample.map <- lapply(samples, function(x){row.names(library.table)[library.table$SAMPLE==x]})

sample.table <- data.frame(row.names=samples, SAMPLE=samples, LIBRARIES=unlist(lapply(sample.map,paste,collapse=",")), LIBCOUNT=unlist(lapply(sample.map,length)))
# Make sure BARCODE (and LINE/SEX/REP if present) is consistent across libraries, pull those into table
for(annot.col in c("BARCODE", intersect(colnames(library.table), c("LINE","SEX","REP")))) {
  if(any(unlist(lapply(sample.map,function(x){length(unique(library.table[x,annot.col]))}))!=1)) {
    cat("Multiple libraries have the same sample ID, but different", annot.col, "values.\n")
    bad.samples <- samples[which(unlist(lapply(sample.map,function(x){length(unique(library.table[x,annot.col]))}))!=1)]
    cat("Problematic samples are:", bad.samples, "\n")
    cat("Further details:\n\n")
    for(bs in bad.samples) {
      print(library.table[library.table$SAMPLE==bs,])
    }
    cat("\n")
    stop("You must correct sample naming conflicts before proceeding!")
    # TO DO: More detailed error message here that identifies the problematice libraries would be useful
  }
  sample.table[,annot.col] <- unlist(lapply(sample.map, function(x){unique(library.table[x,annot.col])}))
}


# If sample_info.txt exists, add that information now
info.file <- "sample_info.txt"
if(file.exists(info.file)) {
  cat("Loading additional sample info from",info.file,"\n")
  info.table <- read.table(info.file, header=T, sep="\t", row.names=1)
  # Warn if there are samples in one table but not the other
  missing.info.rows <- length(setdiff(row.names(sample.table), row.names(info.table)))
  if(missing.info.rows > 0) {
    cat("WARNING: Missing sample info for", missing.info.rows, "samples - these will be NA.\n")
  } else {
    cat("Found info data for all sequenced samples.\n")
  }
  missing.sample.rows <- length(setdiff(row.names(info.table), row.names(sample.table)))
  if(missing.sample.rows > 0) {
    cat("WARNING:",info.file,"contains info for",missing.sample.rows,"that have not been sequenced.\n")
  } else {
    cat("All samples in",info.file,"have been sequenced.\n")
  }
  sample.table <- cbind(sample.table, info.table[row.names(sample.table),])
} else {
  cat(info.file,"is missing, no additional sample info will be added here.\n")
}

# Sum up lib sizes
for(ls.col in c("MICROBE.LIBSZ","REPBASE.LIBSZ","LIBSZ")) {
  sample.table[,ls.col] <- unlist(lapply(sample.map, function(x){sum(library.table[x,ls.col])}))
}

# Add columns for batch structure (designed to be fed directly into EdgeR GLM)
final.batch.cols <- c()
for(b in batches) {
  b.samples <- samples[unlist(lapply(sample.map, function(x){b %in% library.table[x,"BATCH"]}))]
  if(length(final.batch.cols)==1) {
    b.samples.prev <- sample.table[b.samples,final.batch.cols]
  } else {
    b.samples.prev <- apply(sample.table[b.samples,final.batch.cols], 1, sum)
  }
  b.samples.orig <- b.samples[b.samples.prev==0]
  b.samples.reseq <- b.samples[b.samples.prev>0]
  if(length(b.samples.orig) > 0) {
    final.batch.cols <- c(final.batch.cols, b)
    sample.table[,b] <- 0
    sample.table[b.samples.orig,b] <- 1
  }
  if(length(b.samples.reseq) > 0) {
    reseq.col <- paste0(b, "_Reseq")
    final.batch.cols <- c(final.batch.cols, reseq.col)
    sample.table[,reseq.col] <- 0
    sample.table[b.samples.reseq,reseq.col] <- 1
  }
}

cat("Final batch columns and sample counts:\n")
print(apply(sample.table[,final.batch.cols], 2, sum))
cat("\n")

# Progress Summary and Sample Clean-Up

# TO DO: Could add a generalized sample name correction step here, or each project could have its own "correct_sample_table.R" script for that purpose
# But better solution is to correct sample names in the earlier steps

# (OPTIONAL) Load sample log file
# For this part of code to run, create the sample_log_info.txt file with additional info for each sample
# This script will then check to make sure all samples match up, and pull in additional useful columns for downstream technical factor checks
sample.log.file <- "sample_log_info.txt"
if(file.exists(sample.log.file)) {
  cat("Loading additional sample info from",sample.log.file,"and checking against batch-derived info.\n")
  log.table <- read.table(sample.log.file, header=T, sep="\t", row.names=1, na.strings="NA")
  # First identify the sick lines, none of these should be in sample table
  sick.lines <- row.names(log.table)[is.na(log.table$Barcode.Sequence)]
  stopifnot(all(log.table[sick.lines,"NOTES"]=="Sick Line"))
  if(any(sick.lines %in% row.names(sample.table))) {
    stop("Sample table contains one or more sick line samples: ", intersect(sick.lines, row.names(sample.table)))
  }
  log.table <- log.table[setdiff(row.names(log.table),sick.lines),]
  if(!all(row.names(log.table) %in% row.names(sample.table))) {
    cat("WARNING: Have log data for samples without sequencing data:",setdiff(row.names(log.table), row.names(sample.table)))
  }
  if(!all(row.names(sample.table) %in% row.names(log.table))) {
    stop("Missing log data for one or more samples: ", setdiff(row.names(sample.table), row.names(log.table)))
  }

  log.table <- log.table[row.names(sample.table),]
  stopifnot(sum(is.na(log.table))==0)

  # Make sure barcode columns match up between tables
  if(any(sample.table$BARCODE != log.table$Barcode.Sequence)) {
    stop("Some barcodes do not match on log table: ", row.names(sample.table)[sample.table$BARCODE != log.table$Barcode.Sequence])
  }

  # Pull in additional columns of interest
  sample.table[,"PLATE"] <- log.table$Library.Plate
  sample.table[,"WELL"] <- log.table$Library.Well
} else {
  cat(sample.log.file,"does not exist, skipping log file checks.\n")
}

# TO DO:
# DGRP_Baseline_Align version had additional checks making sure each line had two replicates for both sexes
# That code might be useful in another script, but is only needed for larger projects with lots of lines

sample.file <- "sample_master_table.txt"
cat("\nWriting sample-level info table to:",sample.file,"\n")
write.table(sample.table, sample.file, row.names=F, quote=F, sep="\t")

# TO DO: The script for building resequence list should work off this table instead 
# (it will greatly simplify that processing, and there's no point resequencing libraries we aren't going to use!)
