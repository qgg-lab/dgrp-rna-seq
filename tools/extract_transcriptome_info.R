#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 5/16/16
# 
# extract_transcriptome_info.R - Create a table of info for each gene feature in a GFF/GTF file
# Adapted from similar script ~/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/extract_transcriptome_info.R
# 
# TO DO: Consolidate basic functions into a single R script that can both be loaded as a general library
# AND run from the command line using the format: (Rscript.R) (COMMAND) (PARAMS...)
#
# USAGE: extract_transcriptome_info.R transcriptome.gtf
# OUTPUT: transcriptome-gene-info.txt
#


# --- SUBROUTINE --- #

# These functions were copied from  Figaro/Rlib/rnaseq.R on 5/16/16

# Merge function - takes 3 column bed data and does merge operations
mergeExonsGTF <- function(bed.table) {
  # First order by start exon
  bed.table <- bed.table[order(bed.table[,2], decreasing=F),]
  # Then order by chromosome
  bed.table <- bed.table[order(bed.table[,1], decreasing=F),]
  
  # Now loop over each row of bed table up to next to last
  i <- 1
  while(i < nrow(bed.table)) {
    # Check if exon at i and exon at i+1 overlap
    if(bed.table[i,1] == bed.table[i+1,1]) {
      if(bed.table[i,3] >= bed.table[i+1,2]) {
        # If there's an overlap, merge these two into row i, and remove row i+1
        bed.table[i,3] <- max(bed.table[i,3], bed.table[i+1,3])
        keep.rows <- 1:i
        if((i+1) < nrow(bed.table)) {
          keep.rows <- c(keep.rows,(i+2):nrow(bed.table))
        }
        bed.table <- bed.table[keep.rows,]
        # Remind i (gets iterated at end of loop)
        i <- i-1
      }
    }
    i <- i + 1
  }
  
  return(bed.table)
}

# Function to calculate total covered area of merged exons
exonCovgGTF <- function(bed.table) {
  sum(bed.table[,3] - bed.table[,2])
}

# Function to get minimum UCSC span(s) covering all exons
geneSpanGTF <- function(bed.table) {
  geneBody <- c()
  # Loop over distinct chromosomes
  for(chr in unique(bed.table[,1])) {
    start <- min(bed.table[bed.table[,1]==chr,2])
    end <- max(bed.table[bed.table[,1]==chr,3])
    geneBody <- c(geneBody, paste(chr,":",start,"-",end,sep=""))
  }
  return(paste(geneBody, collapse=" // "))
}


# Function to tabulate gene-level feature information from a GTF table
# If using a UCSC GTF file (e.g. refGene), make sure to run the correction function above first
# Then run this function on the corrected version
# Most importantly, this computes the base pairs covered by all exons for each gene symbol, to allow for FPKM normalization
# Set info.file to a file name to auto-save the final table
# UPDATED ON 11/18/14 to handle FlyBase GFF file - haven't confirmed backwards compatibility for RefGene GTF files, but I attempted to maintain this compatibility
# The main change that will be different is now the transcript column is called "TRANSCRIPT" instead of "REFSEQ" for more generality
# UPDATED ON 3/6/17 to also include gene strand!
geneInfoGTF <- function(gtf.file, info.file="") {
  
  # Now load the GTF file
  refgene.gtf <- read.table(gtf.file, sep="\t", stringsAsFactors=F, quote="", comment.char="#")
  # Warning: some of the entries in FlyBase GFF have '#' in the alias but these are mostly transposable_elements and a few gene entries
  # The exon entries seem to be fine, so I'm not worrying about it for now
  # Technically all comment lines start with "##" but there's no way to set comment.char to "##", so it would need to be piped through a grep filter before loading here
  
  # Filter the GTF file down to exon features only
  refgene.exon.gtf <- refgene.gtf[refgene.gtf[,3] == "exon",]
  
  # Now split out the extended key=value pairs in 9th column
  fieldPairs <- strsplit(refgene.exon.gtf[,9], split="[[:space:]]*;[[:space:]]*", fixed=F)
  # Two different formats possible here, one is key=value; the other is key "value";
  fieldPairs <- lapply(fieldPairs, function(x){
    # Trim any leftover leading/trailing white space
    x <- sub("^[[:space:]]*", "", x)
    x <- sub("[[:space:]]*$", "", x)
    pairList=strsplit(x, split="[= ][\"]*", fixed=F)
    newValues=c(lapply(pairList, function(pair){pair[2]}),recursive=T)
    # Trim trailing " from newValues if present
    newValues=sub("\"$", "", newValues)
    names(newValues)=c(lapply(pairList, function(pair){pair[1]}),recursive=T)
    return(newValues)
  })
  
  # Now split out the gene_id and transcript_id
  refgene.exon.gtf[,"GENE"] <- c(lapply(fieldPairs, function(x){x["gene_id"]}),recursive=T)
  refgene.exon.gtf[,"TRANSCRIPT"] <- c(lapply(fieldPairs, function(x){x["transcript_id"]}),recursive=T)
  
  # If all transcripts are NA, then need to use the Parent field instead
  if(all(is.na(refgene.exon.gtf[,"TRANSCRIPT"]))) {
    refgene.exon.gtf[,"TRANSCRIPT"] <- c(lapply(fieldPairs, function(x){x["Parent"]}),recursive=T)
  }
  
  # Now make a list of all unique genes
  refgene.features <- unique(refgene.exon.gtf[,"GENE"])
  # Make sure there are no exons with multiple genes listed
  stopifnot(length(grep(",", refgene.features, fixed=T)) == 0)
  stopifnot(length(grep(" // ", refgene.features, fixed=T)) == 0)
  length(refgene.features)
  
  # Make a table with the following fields:
  feature.info <- data.frame(GENE=refgene.features, stringsAsFactors=F)
  feature.info[,"TRANSCRIPT"] <- as.character(NA)	# A list of corresponding refseq IDs, separated by " // "
  feature.info[,"TRANSCRIPT_COUNT"] <- as.integer(NA)	# Number of distinct transcripts mapping to this gene symbol
  feature.info[,"EXON_COUNT"] <- as.integer(NA)	# Number of total exons in all transcripts mapping to this gene symbol
  feature.info[,"MERGED_EXON_COUNT"] <- as.integer(NA)	# Number of merged exons from all transcripts
  feature.info[,"LENGTH"] <- as.integer(NA)	# Total number of bp covered by all merged exons
  feature.info[,"BODY"] <- as.character(NA)	# Total genomic region covered by all transcripts for this gene
  feature.info[,"STRAND"] <- as.character(NA) # Which strand is gene encoded on? ('.' indicates unstranded or exons on multiple strands)
  feature.info[,"MULTI_CHROM"] <- F	# Flag for when this gene maps to multiple chromosomes, hopefully there will be no such cases
  
  # Now loop over each row and fill in all fields based on the exon gtf data
  # TO DO: This could be drastically sped up with %dopar% / doMC package!
  for(i in 1:nrow(feature.info)) {
    if((i %% 100) == 0) {
      cat(".")
    }
    if((i %% 1000) == 0) {
      cat(i)
    }
    i.gene <- feature.info[i,"GENE"]
    # Extract all exons
    i.exon.gtf <- refgene.exon.gtf[refgene.exon.gtf[,"GENE"]==i.gene,]
    feature.info[i,"EXON_COUNT"] <- nrow(i.exon.gtf)
    # Get info on individual transcripts - make sure to split out multiple transcripts if packed together here
    i.transcript <- unique(i.exon.gtf[,"TRANSCRIPT"])
    i.transcript <- unique(c(strsplit(i.transcript, split=",", fixed=T),recursive=T))
    i.transcript <- unique(c(strsplit(i.transcript, split=" // ", fixed=T),recursive=T))
    feature.info[i,"TRANSCRIPT"] <- paste(i.transcript, collapse=" // ")
    feature.info[i,"TRANSCRIPT_COUNT"] <- length(i.transcript)
    # Now merge all exons
    # First check if exons come from multiple chromosomes
    if(length(unique(i.exon.gtf[,1])) > 1) {
      feature.info[i,"MULTI_CHROM"] <- T
    }
    # Merge sub-routine
    i.exon.merge <- mergeExonsGTF(i.exon.gtf[,c(1,4,5)])
    feature.info[i,"MERGED_EXON_COUNT"] <- nrow(i.exon.merge)
    feature.info[i,"LENGTH"] <- exonCovgGTF(i.exon.merge)
    feature.info[i,"BODY"] <- geneSpanGTF(i.exon.merge)
    i.exon.strands <- unique(i.exon.gtf[,7])
    if(length(i.exon.strands)==1) {
      feature.info[i,"STRAND"] <- i.exon.strands
    } else {
      feature.info[i,"STRAND"] <- "."
    }
  }
  cat("\n")
  
  sum(feature.info[,"MULTI_CHROM"])
  sum(is.na(feature.info))
  
  # Return this info as a table
  if(info.file != "") {
    write.table(feature.info, info.file, row.names=F, col.names=T, sep="\t", quote=F)
  }
  
  return(invisible(feature.info))
  
}


# --- MAIN ROUTINE FOR COMMAND-LINE EXECUTION --- #

# Pull in command line arguments for executing in Rscript mode
# If there are any command-line arguments, parse them and run geneInfoGTF
# If there are no arguments here, script exits silently (so can load the functions using source command)
command.args <- commandArgs(trailingOnly=T)
# TEMP FOR TESTING:
# command.args <- c("transcriptome.gtf")
if(length(command.args) > 0) {
  options(stringsAsFactors=F)

  # Accepts exactly 1 parameter for now
  if(length(command.args) > 1) {
    stop("Too many parameters! This script takes only one command-line parameter: the name of the GFF/GTF file to extract info")
  }
  
  gtf.file <- command.args[1]
  if(!grepl("[.][gG][tTfF][fF]$", gtf.file)) {
    stop("The input file must be a .gtf or .gff file")
  }
  
  if(!file.exists(gtf.file)) {
    stop("Can not find ", gtf.file)
  }
  
  cat("Generating gene info for", gtf.file, "\n")
  file.stub <- sub("[.][gG][tTfF][fF]$", "", gtf.file)
  gene.info.file <- paste0(file.stub, "-gene-info.txt")
  
  geneInfoGTF(gtf.file, info.file=gene.info.file)
  cat("Wrote gene info table to", gene.info.file, "\n")
}
