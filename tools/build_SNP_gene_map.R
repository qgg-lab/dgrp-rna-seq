#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 2/28/17
#
# build_SNP_gene_map.R
#
# Usage:
# eQTL_summary.R EQTL=fdr.results.txt [Options]
#  SNP=     Specifies the list of SNP IDs being used by plink
#  GTF=     path to file of gene models (GTF format)
#           Default: ~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined_gene_models.gtf
#  OUTPUT=  Full path to output file (default to (SNP).gene.map)
#  DIST=    Cutoff for assigning SNPs to genes, in basepairs (Default: 1000)
#  CORES=   Number of threads to use for parallelization
#  
#  
# The purpose of this script is just to build a table linking SNPs to all proximal genes
# Separate row for each SNP,Gene pair, with details of dist and part of gene the SNP falls in
#

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

# For parallelization
library(doMC)

options(stringsAsFactors=F)

# -- Process command-line parameters -- #
usageStr="USAGE:\nbuild_SNP_gene_map.R SNP=freeze2.200line.common.snp [OPTIONS]"
my.args <- commandArgs(trailingOnly=T)
# TEMP TESTING:
# my.args <- c("SNP=freeze2.200line.common.snp", "GTF=known_all_novel_genes/combined_gene_models.gtf", "CORES=4")

if(length(my.args) == 0) {
  stop("Requires at least one argument. ", usageStr)
} else {
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

# SNP=  The file containing line mean expression values to filter and prep for downstream analysis
# Make sure SNP param exists
if(!("SNP" %in% names(my.args))) {
  stop("Missing SNP parameter. ", usageStr)
}

# GTF=  path to file of gene models (GTF format)
# Default: ~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined_gene_models.gtf
setDefault("GTF", "~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined_gene_models.gtf")

# Make sure both files exist before proceeding
if(!file.exists(my.args$SNP)) {
  stop("Missing SNP file: ", my.args$SNP)
}

if(!file.exists(my.args$GTF)) {
  stop("Missing GTF file: ", my.args$GTF)
}

default.output <- sub("[.]snp$", "", my.args$SNP)
default.output <- paste0(default.output, ".snp.gene.map")
setDefault("OUTPUT",default.output)

outDir <- dirname(my.args$OUTPUT)
if(!file.exists(outDir)) {
  dir.create(path=outDir, recursive=T, showWarnings = F)
}

setDefault("DIST",as.integer(1000))

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
  my.args$CORES <- 1
  cat("Setting CORES = 1 by default.\n")
}
# If CORES = ALL, use detectCores()-1
if(is.character(my.args$CORE) & (my.args$CORE == "ALL")) {
  my.args$CORES <- as.integer(detectCores()-1)
  cat("Setting CORES =", my.args$CORES, "\n")
}
# Otherwise, attempt to convert to integer
if(!is.integer(my.args$CORE)) {
  my.args$CORE <- as.integer(my.args$CORE)
  cat("Setting CORES =", my.args$CORE, "\n")
}

# Set up parallel backend with desired number of CPUs
cat("Using", my.args$CORES, "CPUs for parallel tasks.\n")
registerDoMC(my.args$CORES)
cat("\n")


# --- Load Input Files --- #

# Load the SNP File
cat("Loading SNP coordinates from", my.args$SNP, "...\n")
all.snps <- read.table(my.args$SNP)[,1]

# Extract CHR and POS from SNP name
snp.split <- strsplit(all.snps, split="_", fixed=T)
snp.table <- data.frame(
  row.names=all.snps,
  CHR=unlist(lapply(snp.split, "[", 1)),
  POS=as.integer(unlist(lapply(snp.split, "[", 2)))
)
cat("Loaded", nrow(snp.table), "SNP coordinates.\n\n")

# Load the GTF File
cat("Loading gene model data from", my.args$GTF, "\n")
gtf.table <- read.table(my.args$GTF, sep="\t")
colnames(gtf.table) <- c("CHR", "SRC", "TYPE", "START", "END", "SCORE", "STRAND", "FRAME", "ATTR")
gtf.table <- gtf.table[gtf.table$TYPE=="exon",]

# Parse attributes to get gene ID
cat("...parsing attribute fields to get gene IDs...\n")
# TO DO: Use %dopar% to speed this up?
gtf.table.attr <- strsplit(gtf.table$ATTR, split="[;][ ]?")
gtf.table.attr <- lapply(gtf.table.attr, function(x){
  x.split <- strsplit(x, "[ =]", fixed=F)
  x.names <- unlist(lapply(x.split, "[", 1))
  x <- unlist(lapply(x.split, "[", 2))
  x <- sub("\"$", "", sub("^\"", "", x))
  names(x) <- x.names
  return(x)
})
gtf.table.attr.fields <- unique(unlist(lapply(gtf.table.attr, names)))
stopifnot(length(gtf.table.attr.fields) < 20)
gtf.table.attr <- lapply(gtf.table.attr, function(x){x[gtf.table.attr.fields]})
gtf.table.attr <- as.data.frame(matrix(unlist(gtf.table.attr), ncol=length(gtf.table.attr.fields), byrow=T))
colnames(gtf.table.attr) <- gtf.table.attr.fields
stopifnot("gene_id" %in% colnames(gtf.table.attr))
gtf.table <- cbind(gtf.table[,1:8], GENE=gtf.table.attr$gene_id)

drop.exons <- is.na(gtf.table$GENE)
if(any(drop.exons)) {
  cat("Dropping", sum(drop.exons), "exons with no Gene ID.\n")
  gtf.table <- gtf.table[!drop.exons,]
}

cat("Loaded coordinates for", nrow(gtf.table), "exons from", length(unique(gtf.table$GENE)), "genes.\n\n")

# Compute bounds of each gene body
cat("Computing gene body boundaries...\n")
all.genes <- unique(gtf.table$GENE)
# TO DO: Use %dopar% to speed these steps up?
gene.exons <- foreach(gene=all.genes) %dopar% {
  gtf.table[gtf.table$GENE==gene,]
}
names(gene.exons) <- all.genes
gene.chr <- lapply(gene.exons, function(exons){
  unique(exons$CHR)
})
if(any(unlist(lapply(gene.chr, length)) != 1)) {
  stop(sum(unlist(lapply(gene.chr, length)) != 1), "genes have exons on more than one chromosome, cannot properly map gene bodies")
}
gene.strand <- unlist(lapply(gene.exons, function(exons){
  g.strand <- unique(exons$STRAND)
  if(length(g.strand) != 1) {
    return(".")
  } else {
    return(g.strand)
  }
}))
gene.table <- data.frame(
  row.names=all.genes, 
  CHR=unlist(gene.chr),
  START=unlist(lapply(gene.exons, function(exons){min(exons$START)})),
  END=unlist(lapply(gene.exons, function(exons){max(exons$END)})),
  STRAND=gene.strand
)
cat("\n")


# --- Mapping --- #

# First, split all SNPs/Genes by chr (drop those on chr NOT containing SNPs)
snp.chr.tables <- lapply(unique(snp.table$CHR), function(chr){
  snp.table[snp.table$CHR==chr,]
})
names(snp.chr.tables) <- unique(snp.table$CHR)
cat("SNPs cover", length(snp.chr.tables), "chromosomes.\n")

gene.chr.tables <- lapply(unique(gene.table$CHR), function(chr){
  gene.table[gene.table$CHR==chr,]
})
names(gene.chr.tables) <- unique(gene.table$CHR)
cat("Genes cover", length(gene.chr.tables), "chromosomes.\n")

common.chr <- intersect(names(snp.chr.tables), names(gene.chr.tables))
cat(length(common.chr), "chromosomes contain both SNP and genes, these are the only chromosomes that will be processed.\n\n")

# TESTING
# snp.chr.tables <- lapply(snp.chr.tables, head, n = 1000)
# gene.chr.tables <- lapply(gene.chr.tables, head, n = 100)


# Use dopar to loop over GENES, for each gene get list of SNPs w/in gene bounds +/- DIST
snp.gene.map <- foreach(chr=common.chr, .combine='rbind') %do% {
  cat("Mapping SNPs to Genes on", chr, "\n")
  foreach(gene=row.names(gene.chr.tables[[chr]]), .combine='rbind') %dopar% {
    gene.snps <- snp.chr.tables[[chr]]
    gene.snps <- gene.snps[(gene.snps$POS >= (gene.table[gene,"START"]-my.args$DIST)) & (gene.snps$POS <= (gene.table[gene,"END"]+my.args$DIST)),]
    if(nrow(gene.snps) > 0) {
      gene.snps[gene.snps$POS < gene.table[gene,"START"],"DIST"] <- gene.table[gene,"START"] - gene.snps[gene.snps$POS < gene.table[gene,"START"],"POS"]
      gene.snps[(gene.snps$POS >= gene.table[gene,"START"]) & (gene.snps$POS <= gene.table[gene,"END"]),"DIST"] <- 0
      gene.snps[gene.snps$POS > gene.table[gene,"END"],"DIST"] <- gene.snps[gene.snps$POS > gene.table[gene,"END"],"POS"] - gene.table[gene,"END"]
      # Classify those which are 3PRIME, EXON, INTRON, 5PRIME
      gene.snps$PART <- "INTRON"
      if(gene.table[gene,"STRAND"] == "-") {
        gene.snps[gene.snps$POS < gene.table[gene,"START"],"PART"] <- "3PRIME"
        gene.snps[gene.snps$POS > gene.table[gene,"END"],"PART"] <- "5PRIME"
      } else {
        gene.snps[gene.snps$POS < gene.table[gene,"START"],"PART"] <- "5PRIME"
        gene.snps[gene.snps$POS > gene.table[gene,"END"],"PART"] <- "3PRIME"
      }
      my.exons <- gene.exons[[gene]]
      overlap.exons <- unlist(lapply(gene.snps$POS, function(x){
        any((my.exons$START <= x) & (my.exons$END >= x))
      }))
      gene.snps[overlap.exons,"PART"] <- "EXON"
      return(data.frame(
        SNP.CHR=gene.snps$CHR,
        SNP.POS=gene.snps$POS,
        SNP=row.names(gene.snps),
        GENE=rep(gene, times=nrow(gene.snps)),
        DIST=gene.snps$DIST,
        GENE.PART=gene.snps$PART,
        GENE.STRAND=rep(gene.table[gene,"STRAND"], times=nrow(gene.snps))
      ))
    } else {
      return(data.frame(SNP.CHR=character(), SNP.POS=integer(), SNP=character(), GENE=character(), DIST=integer(), GENE.PART=character(), GENE.STRAND=character()))
    }
  }
}

# Sort by SNP.CHR, SNP.POS, then drop these columns
cat("Sorting SNP-GENE connections by SNP coordinate...\n")
snp.gene.map <- snp.gene.map[order(snp.gene.map$SNP.POS),]
snp.gene.map <- snp.gene.map[order(snp.gene.map$SNP.CHR),]
snp.gene.map <- snp.gene.map[,3:ncol(snp.gene.map)]
cat("Finished processinging all SNPs and Genes.\n\n")

# Output table
cat("Writing SNP-Gene map to:", my.args$OUTPUT, "\n")
write.table(snp.gene.map, my.args$OUTPUT, row.names=F, sep="\t", quote=F)
cat("\n")

cat("Script completed successfully!\n\n")
print(proc.time())
cat("\n")
