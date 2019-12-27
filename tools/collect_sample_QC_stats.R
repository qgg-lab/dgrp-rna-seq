#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 12/12/14
#
# Collect sample QC stats into a single table
#
# USAGE:
# collect_sample_QC_stats.R BATCH1 [BATCH2 ...]
#

# setwd("~/Projects/DGRP_Baseline_RNAseq")

library(matrixStats)

options(stringsAsFactors=F)

# Take the list of batches to analyze from command line params
usageStr <- "USAGE:\nRscript collect_sample_QC_stats.R BATCH1 [BATCH2 ...]"
myArgs <- commandArgs(trailingOnly=T)
if(length(myArgs) == 0) {
  # Can't run without at least one batch!
  stop("Must specify one or more batches to tabulate QC stats. ", usageStr)
} else {
  all.batches <- myArgs
}

# Remove trailing .txt from these files (just in case)
all.batches <- sub("[.]txt$", "", all.batches)

# For interactive code testing:
# all.batches <- c("batch_141031","batch_150318","batch_150427","batch_150527","batch_150707","batch_150917")


# Loop over the code below for each batch, then rbind the tables together at the end
complete.table <- NULL
for(batch in all.batches) {
  cat("Processing",batch,"\n")
  batch.file <- paste(batch, ".txt", sep="")
  
  batch.table <- read.table(batch.file, sep="\t")
  if(ncol(batch.table) == 7) {
  	colnames(batch.table) <- c("SAMPLE.ID","SAMPLE.NAME","FLOWCELL","LANE","INDEX","BARCODE","SEQ.DATE")
  } else {
    colnames(batch.table) <- c("SAMPLE.ID","SAMPLE.NAME","FLOWCELL","LANE","BARCODE","SEQ.DATE")
  }
  row.names(batch.table) <- batch.table[,"SAMPLE.ID"]
  batch.table[,"BATCH"] <- batch
  # NOTE: Original DGRP Baseline version of this script also parsed out line #, sex, and replicate # here
  # But this format varies by project, and is not actually used in this script, so it was taken out
  
  # Now start collecting sample information
  
  # Starting library size (ideally would like to know even before any QC filtering done by Casava, but the best I can do is read count from the fastq file I started with)
  # Use results from initial fastqc runs
  batch.table[,"TOTAL.READS"] <- as.integer(NA)
  for(sampleID in row.names(batch.table)) {
    fastqc.data.file <- paste("fastqc/",batch,"/",sampleID,"_fastqc/fastqc_data.txt",sep="")
    if(!file.exists(fastqc.data.file)) {
      stop("MISSING: ",fastqc.data.file)
    }
    cmd <- paste("grep '^Total Sequences'",fastqc.data.file,"| awk '{print $3}'")
    cmd.pipe <- pipe(description=cmd, open="r")
    cmd.result <- scan(cmd.pipe, quiet=T)
    close(cmd.pipe)
    stopifnot(length(cmd.result) == 1)
    batch.table[sampleID,"TOTAL.READS"] <- as.integer(cmd.result)
  }
  stopifnot(sum(is.na(batch.table[,"TOTAL.READS"])) == 0)
  
  # Starting library encoding
  batch.table[,"ENCODING"] <- as.character(NA)
  for(sampleID in row.names(batch.table)) {
    fastqc.data.file <- paste("fastqc/",batch,"/",sampleID,"_fastqc/fastqc_data.txt",sep="")
    if(!file.exists(fastqc.data.file)) {
      stop("MISSING: ",fastqc.data.file)
    }
    cmd <- paste("grep '^Encoding'",fastqc.data.file,"| cut -f 2")
    cmd.pipe <- pipe(description=cmd, open="r")
    cmd.result <- readLines(cmd.pipe)
    close(cmd.pipe)
    stopifnot(length(cmd.result) == 1)
    batch.table[sampleID,"ENCODING"] <- cmd.result
  }
  stopifnot(sum(is.na(batch.table[,"ENCODING"])) == 0)
  
  # Percent of reads removed by cutadapt
  # Pull the "ShortReads" column out of trimmed/batch_141031.cutadapt.summary.txt and divide that into total reads
  cutadapt.summary.file <- paste("trimmed/", batch, ".cutadapt.summary.txt", sep="")
  stopifnot(file.exists(cutadapt.summary.file))
  cutadapt.summary.table <- read.table(cutadapt.summary.file, header=T, sep="\t", comment.char="", row.names=1)
  row.names(cutadapt.summary.table) <- sub("_trimmed.fastq.gz$", "", row.names(cutadapt.summary.table))
  stopifnot(all(row.names(batch.table) %in% row.names(cutadapt.summary.table)))
  cutadapt.summary.table <- cutadapt.summary.table[row.names(batch.table),]
  batch.table[,"FRAC.READS.TRIMMED"] <- cutadapt.summary.table[,"ShortReads"] / batch.table[,"TOTAL.READS"]
  
  # Percent of bases removed by cutadapt
  # This will have to be pulled from cutadapt log files which have uninformative file names (doesn't say which sample)
  batch.table[,"FRAC.BASES.TRIMMED"] <- as.numeric(NA)
  for(log.file in dir(paste("trimmed/",batch,sep=""), pattern="^run_cutadapt_slurm_.*log$", full.names=T)) {
    # Load the complete log file
    log.lines <- readLines(log.file)
    # Figure out which sample this corresponds to
    log.run.line <- grep('^Running:', log.lines, value=T)
    stopifnot(length(log.run.line) == 1)
    log.run.tokens <- unlist(strsplit(log.run.line, split="[[:space:]]+"))
    log.out.file <- log.run.tokens[which(log.run.tokens=="-o")[1]+1]
    sampleID <- sub("_trimmed.fastq.gz","",sub(paste("^.*/",batch,"/",sep=""),"",log.out.file))
    stopifnot(sampleID %in% row.names(batch.table))
    stopifnot(is.na(batch.table[sampleID,"FRAC.BASES.TRIMMED"]))
    
    # Now figure out the % of bases trimmed
    log.base.line <- grep('Trimmed bases:', log.lines, value=T)
    stopifnot(length(log.base.line) == 1)
    log.base.tokens <- unlist(strsplit(log.base.line, split="[[:space:]]+"))
    log.base.perc.str <- grep("^\\([0-9.]+%$", log.base.tokens, value=T)
    stopifnot(length(log.base.perc.str) == 1)
    log.base.perc <- sub("^\\(","",sub("%$","",log.base.perc.str))
    batch.table[sampleID,"FRAC.BASES.TRIMMED"] <- as.numeric(log.base.perc) / 100
  }
  stopifnot(sum(is.na(batch.table[,"FRAC.BASES.TRIMMED"])) == 0)
  
  # Percent of reads removed by ribo-filter
  # Load the read counts for trimmed and ribo-filtered fastq files
  trimmed.fastq.wcl <- read.table("trimmed/fastq.gz.stats", sep="\t")
  colnames(trimmed.fastq.wcl) <- c("FILE","LINE.COUNT","FILE.SIZE","MD5.SUM")
  trimmed.fastq.wcl[,"READ.COUNT"] <- as.integer(trimmed.fastq.wcl[,"LINE.COUNT"] / 4)
  trimmed.fastq.wcl <- trimmed.fastq.wcl[grep(batch,trimmed.fastq.wcl[,"FILE"]),]
  row.names(trimmed.fastq.wcl) <- sub(paste("^.*/",batch,"/",sep=""),"",sub("_trimmed.fastq.gz","",trimmed.fastq.wcl[,"FILE"]))
  stopifnot(all(row.names(batch.table) %in% row.names(trimmed.fastq.wcl)))
  trimmed.fastq.wcl <- trimmed.fastq.wcl[row.names(batch.table),]
  
  ribofiltered.fastq.wcl <- read.table("ribofiltered/fastq.gz.stats", sep="\t")
  colnames(ribofiltered.fastq.wcl) <- c("FILE","LINE.COUNT","FILE.SIZE","MD5.SUM")
  ribofiltered.fastq.wcl[,"READ.COUNT"] <- as.integer(ribofiltered.fastq.wcl[,"LINE.COUNT"] / 4)
  ribofiltered.fastq.wcl <- ribofiltered.fastq.wcl[grep(batch,ribofiltered.fastq.wcl[,"FILE"]),]
  row.names(ribofiltered.fastq.wcl) <- sub(paste("^.*/",batch,"/",sep=""),"",sub("_filtered.fastq.gz","",ribofiltered.fastq.wcl[,"FILE"]))
  stopifnot(all(row.names(batch.table) %in% row.names(ribofiltered.fastq.wcl)))
  ribofiltered.fastq.wcl <- ribofiltered.fastq.wcl[row.names(batch.table),]
  
  batch.table[,"FRAC.READS.RIBO"] <- (trimmed.fastq.wcl[,"READ.COUNT"] - ribofiltered.fastq.wcl[,"READ.COUNT"]) / trimmed.fastq.wcl[,"READ.COUNT"]
  stopifnot(sum(is.na(batch.table[,"FRAC.READS.RIBO"])) == 0)
  
  # Percent of reads removed by Microbe-filter
  microbefiltered.fastq.wcl <- read.table("microbefiltered/fastq.gz.stats", sep="\t")
  colnames(microbefiltered.fastq.wcl) <- c("FILE","LINE.COUNT","FILE.SIZE","MD5.SUM")
  microbefiltered.fastq.wcl[,"READ.COUNT"] <- as.integer(microbefiltered.fastq.wcl[,"LINE.COUNT"] / 4)
  microbefiltered.fastq.wcl <- microbefiltered.fastq.wcl[grep(batch,microbefiltered.fastq.wcl[,"FILE"]),]
  row.names(microbefiltered.fastq.wcl) <- sub(paste("^.*/",batch,"/",sep=""),"",sub("_filtered.fastq.gz","",microbefiltered.fastq.wcl[,"FILE"]))
  stopifnot(all(row.names(batch.table) %in% row.names(microbefiltered.fastq.wcl)))
  microbefiltered.fastq.wcl <- microbefiltered.fastq.wcl[row.names(batch.table),]
  
  batch.table[,"FRAC.READS.MICROBE"] <- (ribofiltered.fastq.wcl[,"READ.COUNT"] - microbefiltered.fastq.wcl[,"READ.COUNT"]) / ribofiltered.fastq.wcl[,"READ.COUNT"]
  stopifnot(sum(is.na(batch.table[,"FRAC.READS.MICROBE"]))==0)
  
  # Percent of reads removed by RepBase-filter
  repfiltered.fastq.wcl <- read.table("repfiltered/fastq.gz.stats", sep="\t")
  colnames(repfiltered.fastq.wcl) <- c("FILE","LINE.COUNT","FILE.SIZE","MD5.SUM")
  repfiltered.fastq.wcl[,"READ.COUNT"] <- as.integer(repfiltered.fastq.wcl[,"LINE.COUNT"] / 4)
  repfiltered.fastq.wcl <- repfiltered.fastq.wcl[grep(batch,repfiltered.fastq.wcl[,"FILE"]),]
  row.names(repfiltered.fastq.wcl) <- sub(paste("^.*/",batch,"/",sep=""),"",sub("_filtered.fastq.gz","",repfiltered.fastq.wcl[,"FILE"]))
  stopifnot(all(row.names(batch.table) %in% row.names(repfiltered.fastq.wcl)))
  repfiltered.fastq.wcl <- repfiltered.fastq.wcl[row.names(batch.table),]
  
  batch.table[,"FRAC.READS.REPEAT"] <- (microbefiltered.fastq.wcl[,"READ.COUNT"] - repfiltered.fastq.wcl[,"READ.COUNT"]) / microbefiltered.fastq.wcl[,"READ.COUNT"]
  batch.table[,"FILTERED.READS"] <- repfiltered.fastq.wcl[,"READ.COUNT"]
  stopifnot(sum(is.na(batch.table[,c("FRAC.READS.MICROBE","FILTERED.READS")])) == 0)
  
  # Read Length Statistics after filtering - how best to summarize these?  Median?  Mean?  StDev?
  batch.table[,"READ.LEN.MEDIAN"] <- as.numeric(NA)
  batch.table[,"READ.LEN.MEAN"] <- as.numeric(NA)
  batch.table[,"READ.LEN.STDEV"] <- as.numeric(NA)
  batch.table[,"READ.LEN.MODE"] <- as.integer(NA)
  for(sampleID in row.names(batch.table)) {
    fastqc.data.file <- paste("fastqc_repfiltered/",batch,"/",sampleID,"_filtered_fastqc/fastqc_data.txt",sep="")
    stopifnot(file.exists(fastqc.data.file))
    fastqc.data <- readLines(fastqc.data.file)
    length.data.start <- grep("^>>Sequence Length Distribution",fastqc.data)[1]+2
    length.data.end <- grep("^>>Sequence Duplication Levels",fastqc.data)[1]-2
    length.data <- strsplit(fastqc.data[length.data.start:length.data.end], split="\t", fixed=T)
    length.values <- sapply(length.data, function(x){
      if(length(grep("-",x[1])) == 1) {
        range.vals <- as.numeric(unlist(strsplit(x[1], "-")))
        return(mean(range.vals))
      } else {
        return(mean(as.numeric(x[1])))
      }
    })
    length.weights <- as.numeric(sapply(length.data, "[", 2))
    batch.table[sampleID,"READ.LEN.MEDIAN"] <- weightedMedian(x=length.values, w=length.weights)
    batch.table[sampleID,"READ.LEN.MEAN"] <- weighted.mean(x=length.values, w=length.weights)
    batch.table[sampleID,"READ.LEN.STDEV"] <- sqrt(weightedVar(x=length.values, w=length.weights))
    batch.table[sampleID,"READ.LEN.MODE"] <- as.integer(ceiling(length.values[which.max(length.weights)]))
  }
  stopifnot(sum(is.na(batch.table[,c("READ.LEN.MEDIAN","READ.LEN.MEAN","READ.LEN.STDEV")]))==0)
  
  # STAR Alignment Stats:
  # Uniquely Aligned Read Count, 
  # % of Reads Uniquely Aligned,
  # Multiple Aligned Read Count (Doesn't count reads with "too many" alignments),
  # % of Reads w/ Multiple Alignments
  # % of Reads too short
  # NOTE: These %s won't add to 100%, the remaining fraction should be quite small and consists of:
  # reads with too many alignments, reads with too many mismatches, and reads unaligned for "other" reasons
  batch.table[,"UNIQUE.ALIGN.READS"] <- as.integer(NA)
  batch.table[,"FRAC.UNIQUE.ALIGN"] <- as.numeric(NA)
  batch.table[,"MULTI.ALIGN.READS"] <- as.integer(NA)
  batch.table[,"FRAC.MULTI.ALIGN"] <- as.numeric(NA)
  batch.table[,"FRAC.TOO.SHORT"] <- as.numeric(NA)
  for(sampleID in row.names(batch.table)) {
    star.log.file <- paste("star/",batch,"/",sampleID,"/Log.final.out",sep="")
    stopifnot(file.exists(star.log.file))
    extract.cmd <- paste("grep 'Uniquely mapped reads number'",star.log.file,"| awk -F'\t' '{print $2}'")
    extract.cmd <- pipe(description=extract.cmd, open="r")
    extract.val <- scan(extract.cmd, quiet=T)
    close(extract.cmd)
    stopifnot(length(extract.val) == 1)
    batch.table[sampleID,"UNIQUE.ALIGN.READS"] <- as.integer(extract.val)
    batch.table[sampleID,"FRAC.UNIQUE.ALIGN"] <- batch.table[sampleID,"UNIQUE.ALIGN.READS"] / batch.table[sampleID,"FILTERED.READS"]
    
    extract.cmd <- paste("grep 'Number of reads mapped to multiple loci'",star.log.file,"| awk -F'\t' '{print $2}'")
    extract.cmd <- pipe(description=extract.cmd, open="r")
    extract.val <- scan(extract.cmd, quiet=T)
    close(extract.cmd)
    stopifnot(length(extract.val) == 1)
    batch.table[sampleID,"MULTI.ALIGN.READS"] <- as.integer(extract.val)
    batch.table[sampleID,"FRAC.MULTI.ALIGN"] <- batch.table[sampleID,"MULTI.ALIGN.READS"] / batch.table[sampleID,"FILTERED.READS"]
    
    extract.cmd <- paste("grep '% of reads unmapped: too short'",star.log.file,"| awk -F'\t' '{print $2}' | sed 's/%$//'")
    extract.cmd <- pipe(description=extract.cmd, open="r")
    extract.val <- scan(extract.cmd, quiet=T)
    close(extract.cmd)
    stopifnot(length(extract.val) == 1)
    batch.table[sampleID,"FRAC.TOO.SHORT"] <- extract.val / 100
  }
  stopifnot(sum(is.na(batch.table[,c("UNIQUE.ALIGN.READS","FRAC.UNIQUE.ALIGN","MULTI.ALIGN.READS","FRAC.MULTI.ALIGN","FRAC.TOO.SHORT")])) == 0)
  
  # % of Uniquely Aligned Reads in Known Genes
  batch.table[,"KNOWN.GENE.READS"] <- as.integer(NA)
  batch.table[,"FRAC.KNOWN.GENES"] <- as.numeric(NA)
  for(sampleID in row.names(batch.table)) {
    count.file <- paste("htseq/",batch,"/",sampleID,"_STAR_counts.txt",sep="")
    count.table <- read.table(count.file, sep="\t")
    batch.table[sampleID,"KNOWN.GENE.READS"] <- sum(count.table[grep("^__", count.table[,1], invert=T),2])
    batch.table[sampleID,"FRAC.KNOWN.GENES"] <- batch.table[sampleID,"KNOWN.GENE.READS"] / batch.table[sampleID,"UNIQUE.ALIGN.READS"]
  }
  stopifnot(sum(is.na(batch.table[,c("KNOWN.GENE.READS","FRAC.KNOWN.GENES")])) == 0)
  
  complete.table <- rbind(complete.table, batch.table)
}

write.table(complete.table, "sample_QC_stats.txt", row.names=F, sep="\t", quote=F)

cat("Finished successfully!")

proc.time()
