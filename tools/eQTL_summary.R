#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 2/28/17
#
# eQTL_summary.R
#
# Usage:
# eQTL_summary.R EQTL=fdr.results.txt [Options]
#  EQTL=    Specifies the main input file of EQTLs from eQTL_perm_fdr.R
#  MAP=     Specifies the SNP-Gene map file (created by build_SNP_gene_map.R)
#  GENES=   Specifies the Gene info file
#           Default: ~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined-gene-info.txt
#  [DEPRECATED] GTF=     path to file of gene models (GTF format)
#  [DEPRECATED]          Default: ~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined_gene_models.gtf
#  OUTPUT=  Path and file stub for all output files (defaults to EQTL minus the .results.txt suffix)
#  DIST=    Cutoff for assigning SNPs to genes, in basepairs (Default: 1000)
#  CORES=   Number of threads to use for parallelization
#  
#  
# Summarize eQTL results:
# 1) Split eQTL classifications into cis and trans (output in seperate tables?)
#    Classify as cis if: A) SNP w/in 1kb of gene; B) closest gene to SNP is the target gene? (Up to some additional MAX?)
# 2) For trans eQTLs only, map each QTL to closest/best gene, output as list of network edges
#    Assignment rule should be: In gene (how to handle multiple?); 5' Upstream (1kb); 3' Downstream (1kb)
#    Should also consider if any eQTL is a cis eQTL for some other gene?
#


# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

# For parallelization
library(doMC)

options(stringsAsFactors=F)

# -- Process command-line parameters -- #
usageStr="USAGE:\neQTL_summary.R EQTL=fdr.results.txt MAP=snp.gene.map [Options]"
my.args <- commandArgs(trailingOnly=T)
# TEMP TESTING:
# my.args <- c("EQTL=known_all_novel_genes/plink/ExprF.0.05.fdr.results.txt", "MAP=freeze2.200line.common.snp.gene.map", "GENES=known_all_novel_genes/combined-gene-info.txt", "CORES=4")
# DEPRECATED: my.args <- c("EQTL=known_all_novel_genes/plink/ExprF.0.05.fdr.results.txt", "MAP=freeze2.200line.common.snp.gene.map", "GENES=known_all_novel_genes/combined-gene-info.txt", "GTF=known_all_novel_genes/combined_gene_models.gtf", "CORES=4")

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

# EQTL=  The file containing line mean expression values to filter and prep for downstream analysis
# Make sure EQTL param exists
if(!("EQTL" %in% names(my.args))) {
  stop("Missing EQTL parameter. ", usageStr)
}

# MAP=  Specifies the SNP-Gene map file (created by build_SNP_gene_map.R)
# Make sure MAP param is specified
if(!("MAP" %in% names(my.args))) {
  stop("Missing MAP parameter. ", usageStr)
}

# GENES=  Specifies the Gene info file
#         Default: ~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined-gene-info.txt
setDefault("GENES", "~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined-gene-info.txt")

# GTF=  path to file of gene models (GTF format)
# Default: ~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined_gene_models.gtf
# DEPRECATED
# setDefault("GTF", "~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/combined_gene_models.gtf")

# Make sure both files exist before proceeding
if(!file.exists(my.args$EQTL)) {
  stop("Missing EQTL file: ", my.args$EQTL)
}

if(!file.exists(my.args$MAP)) {
  stop("Missing MAP file: ", my.args$MAP)
}

if(!file.exists(my.args$GENES)) {
  stop("Missing GENES file: ", my.args$GENES)
}

# DEPRECATED
# if(!file.exists(my.args$GTF)) {
#  stop("Missing GTF file: ", my.args$GTF)
# }

default.outstub <- sub("[.]txt$", "", my.args$EQTL)
default.outstub <- sub("[.]results$", "", default.outstub)
setDefault("OUTPUT",default.outstub)

outDir <- dirname(my.args$OUTPUT)
if(!file.exists(outDir)) {
  dir.create(path=outDir, recursive=T, showWarnings = F)
}

# DEPRECATED
# setDefault("DIST",as.integer(1000))

# TO DO: REMOVE MULTI-THREADING FOR THIS SCRIPT?
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

# Load the EQTL File
eqtl.table <- read.table(my.args$EQTL, header=T, sep="\t", row.names=1)
cat("Loaded eQTL results for", nrow(eqtl.table), "features with FDR =", max(eqtl.table$FDR), "from", my.args$EQTL, "\n")
# Drop all features with no eQTLs
drop.features <- is.na(eqtl.table$SIGSNPS)
if(any(drop.features)) {
  cat("Dropping", sum(drop.features), "features with no eQTLs.\n")
  eqtl.table <- eqtl.table[!drop.features,]
  cat(nrow(eqtl.table), "features remain for analysis.\n")
}
cat("\n")

# TEMP TESTING: Drop down to 100 Genes
# eqtl.table <- eqtl.table[1:100,]

# Load the MAP File
snp.gene.map <- read.table(my.args$MAP, header=T, sep="\t")
cat("Loaded", nrow(snp.gene.map), "proximity associations between", length(unique(snp.gene.map$SNP)), "SNPs and", length(unique(snp.gene.map$GENE)), "genes.\n")

# Load the GENES File
gene.table <- read.table(my.args$GENES, header=T, sep="\t", row.names=1)
cat("Loaded gene coordinates for", nrow(gene.table), "genes from:", my.args$GENES, "\n")

# DEPRECATED
# # Load the GTF File
# cat("Loading gene model data from", my.args$GTF, "\n")
# gtf.table <- read.table(my.args$GTF, sep="\t")
# colnames(gtf.table) <- c("CHR", "SRC", "TYPE", "START", "END", "SCORE", "STRAND", "FRAME", "ATTR")
# gtf.table <- gtf.table[gtf.table$TYPE=="exon",]
# 
# # Parse attributes to get gene ID
# cat("...parsing attribute fields to get gene IDs...\n")
# # TO DO: Use %dopar% to speed this up?
# gtf.table.attr <- strsplit(gtf.table$ATTR, split="[;][ ]?")
# gtf.table.attr <- lapply(gtf.table.attr, function(x){
#   x.split <- strsplit(x, "[ =]", fixed=F)
#   x.names <- unlist(lapply(x.split, "[", 1))
#   x <- unlist(lapply(x.split, "[", 2))
#   x <- sub("\"$", "", sub("^\"", "", x))
#   names(x) <- x.names
#   return(x)
# })
# gtf.table.attr.fields <- unique(unlist(lapply(gtf.table.attr, names)))
# stopifnot(length(gtf.table.attr.fields) < 20)
# gtf.table.attr <- lapply(gtf.table.attr, function(x){x[gtf.table.attr.fields]})
# gtf.table.attr <- as.data.frame(matrix(unlist(gtf.table.attr), ncol=length(gtf.table.attr.fields), byrow=T))
# colnames(gtf.table.attr) <- gtf.table.attr.fields
# stopifnot("gene_id" %in% colnames(gtf.table.attr))
# gtf.table <- cbind(gtf.table[,1:8], GENE=gtf.table.attr$gene_id)
# 
# drop.exons <- is.na(gtf.table$GENE)
# if(any(drop.exons)) {
#   cat("Dropping", sum(drop.exons), "exons with no Gene ID.\n")
#   gtf.table <- gtf.table[!drop.exons,]
# }
# 
# cat("Loaded coordinates for", nrow(gtf.table), "exons from", length(unique(gtf.table$GENE)), "genes.\n\n")

# Make sure all genes in SNP-Gene map are in gene.table
missing.genes <- setdiff(unique(snp.gene.map$GENE), row.names(gene.table))
if(length(missing.genes) > 0) {
  stop(my.args$GENE, " missing gene coordinates for ", length(missing.genes), " in ", my.args$MAP, " including: ", paste(head(missing.genes), collapse=","))
}

# DEPRECATED
# # Compute bounds of each gene body
# cat("Computing gene body boundaries...\n")
# all.genes <- unique(gtf.table$GENE)
# # TO DO: Use %dopar% to speed these steps up?
# gene.exons <- lapply(all.genes, function(gene){
#   gtf.table[gtf.table$GENE==gene,]
# })
# names(gene.exons) <- all.genes
# gene.chr <- lapply(gene.exons, function(exons){
#   unique(exons$CHR)
# })
# if(any(unlist(lapply(gene.chr, length)) != 1)) {
#   stop(sum(unlist(lapply(gene.chr, length)) != 1), "genes have exons on more than one chromosome, cannot properly map gene bodies")
# }
# gene.strand <- unlist(lapply(gene.exons, function(exons){
#   g.strand <- unique(exons$STRAND)
#   if(length(g.strand) != 1) {
#     return(".")
#   } else {
#     return(g.strand)
#   }
# }))
# gene.bodies <- data.frame(
#   row.names=all.genes, 
#   CHR=unlist(gene.chr),
#   START=unlist(lapply(gene.exons, function(exons){min(exons$START)})),
#   END=unlist(lapply(gene.exons, function(exons){max(exons$END)})),
#   STRAND=gene.strand
# )
# cat("\n")

# Now expand eQTL table so there is a row for every Gene-SNP connection
cat("Expanding eQTL result table into individual Gene-SNP connections...\n")
start.time <- proc.time()
eqtl.snps <- strsplit(eqtl.table$SIGSNPS, split=",", fixed=T)
names(eqtl.snps) <- row.names(eqtl.table)
eqtl.pvals <- strsplit(eqtl.table$SIGPVALS, split=",", fixed=T)
eqtl.pvals <- lapply(eqtl.pvals, as.numeric)
names(eqtl.pvals) <- row.names(eqtl.table)
stopifnot(unlist(lapply(eqtl.snps, length))==unlist(lapply(eqtl.pvals, length)))
eqtl.snp.table <- foreach(gene=row.names(eqtl.table), .combine=rbind) %dopar% {
  data.frame(
    GENE=rep(gene, times=length(eqtl.snps[[gene]])),
    SNP=eqtl.snps[[gene]],
    PVAL=eqtl.pvals[[gene]]
  )
}
end.time <- proc.time()
cat("Completed in:\n")
print(end.time-start.time)
cat("\n")

# Check if all gene IDs in eQTL table are in the Gene table
# Don't cause an error, but this tells us if the features are genes or something else
feature.in.map <- unique(eqtl.snp.table$GENE) %in% row.names(gene.table)
names(feature.in.map) <- unique(eqtl.snp.table$GENE)
if(all(feature.in.map)) {
  cat("All features have gene coordinates, presuming this is gene expression eQTL results.\n")
  gene.mode <- T
} else if(all(!feature.in.map)) {
  cat("No eQTL features have gene coordinates, presuming this is microbe or transposon eQTL results.\n")
  gene.mode <- F
} else {
  cat("WARNING: Missing gene coordinates for", sum(!feature.in.map), "genes in eQTL table, including:", paste(head(names(feature.in.map)[!feature.in.map]), collapse=","), "\n")
  cat("Reporting of gene-SNP distances will be suppressed.\n")
  gene.mode <- F
}

if(gene.mode) {
  # Extract gene coordinates from gene.table
  cat("Mapping gene coordinates onto eQTL table...\n")
  gene.split <- strsplit(gene.table$BODY, split="[:-]")
  gene.table$CHR <- unlist(lapply(gene.split, "[", 1))
  gene.table$START <- as.integer(unlist(lapply(gene.split, "[", 2)))
  gene.table$END <- as.integer(unlist(lapply(gene.split, "[", 3)))
  gene.bodies <- gene.table[,c("CHR","START","END","STRAND")]
  
  # Map Gene Coordinates onto table (When features = Genes only)
  eqtl.snp.table <- cbind(eqtl.snp.table, gene.bodies[eqtl.snp.table$GENE,])
  colnames(eqtl.snp.table)[4:7] <- paste0("GENE.", colnames(eqtl.snp.table)[4:7])
  cat("\n")
}

# Make a list of features that have NO linked SNPs
# And therefore can NOT have cis eQTLs
snpless.features <- setdiff(row.names(eqtl.table), unique(snp.gene.map$GENE))
cat(length(snpless.features), "features have no SNPs mapped for them, these features cannot have cis eQTL.\n\n")

# Drop SNP-Gene associations for all SNPs not in eqtl.snp.table?
map.eqtl <- snp.gene.map$SNP %in% unique(eqtl.snp.table$SNP)
cat("Dropping", sum(!map.eqtl), "SNP-Gene associations for non-eQTL SNPs.\n")
snp.gene.map <- snp.gene.map[map.eqtl,]
cat(nrow(snp.gene.map), "SNP-Gene associations remain for mapping of eQTLs.\n\n")

# Extract CHR and POS from SNP name
# (Always do this, it might be useful for hotspot scanning later, even for microbiome/transposon)
cat("Extracting SNP coordinates...\n")
eqtl.snp.split <- strsplit(eqtl.snp.table$SNP, split="_", fixed=T)
eqtl.snp.table$SNP.CHR <- unlist(lapply(eqtl.snp.split, "[", 1))
eqtl.snp.table$SNP.POS <- as.integer(unlist(lapply(eqtl.snp.split, "[", 2)))
cat("\n")

# Compute TARGET.DIST between gene feature and SNP (gene-mode only)
if(gene.mode) {
  cat("Computing target distance between each eQTL and target gene...\n")
  start.time <- proc.time()
  overlap.rows <- (eqtl.snp.table$GENE.CHR == eqtl.snp.table$SNP.CHR) & (eqtl.snp.table$SNP.POS >= eqtl.snp.table$GENE.START) & (eqtl.snp.table$SNP.POS <= eqtl.snp.table$GENE.END)
  dist.rows <- (eqtl.snp.table$GENE.CHR == eqtl.snp.table$SNP.CHR) & !overlap.rows
  diff.chr.rows <- eqtl.snp.table$GENE.CHR != eqtl.snp.table$SNP.CHR
  eqtl.snp.table[,"TARGET.DIST"] <- as.integer(NA)
  # TO DO: Could speed this up further with dopar?
  eqtl.snp.table[dist.rows,"TARGET.DIST"] <- apply(cbind(
    abs(eqtl.snp.table[dist.rows,"GENE.START"]-eqtl.snp.table[dist.rows,"SNP.POS"]),
    abs(eqtl.snp.table[dist.rows,"GENE.END"]-eqtl.snp.table[dist.rows,"SNP.POS"])
  ), 1, min)
  eqtl.snp.table[diff.chr.rows,"TARGET.DIST"] <- Inf
  eqtl.snp.table[overlap.rows,"TARGET.DIST"] <- 0
  end.time <- proc.time()
  stopifnot(sum(is.na(eqtl.snp.table$TARGET.DIST))==0)
  cat("Finished in:\n")
  print(end.time-start.time)
  cat("\n")
  
  # DEPRECATED:
  # eqtl.snp.table$TARGET.DIST <- foreach(i=1:nrow(eqtl.snp.table), .combine='c') %dopar% {
  #   if(eqtl.snp.table[i,"GENE.CHR"] == eqtl.snp.table[i,"SNP.CHR"]) {
  #     if((eqtl.snp.table[i,"SNP.POS"] >= eqtl.snp.table[i,"GENE.START"]) & (eqtl.snp.table[i,"SNP.POS"] <= eqtl.snp.table[i,"GENE.END"])) {
  #       0
  #     } else {
  #       min(abs(eqtl.snp.table[i,"GENE.START"]-eqtl.snp.table[i,"SNP.POS"]), abs(eqtl.snp.table[i,"GENE.END"]-eqtl.snp.table[i,"SNP.POS"]))
  #     }
  #   } else {
  #     Inf
  #   }
  # }
  # cat("\n")
}

# Use Gene.SNP as row name in both eqtl.snp.table and gene.snp.map
row.names(eqtl.snp.table) <- paste(eqtl.snp.table$GENE, eqtl.snp.table$SNP, sep=".")
row.names(snp.gene.map) <- paste(snp.gene.map$GENE, snp.gene.map$SNP, sep=".")


# Separate cis and trans eQTLs
# Cis eQTLs are strictly defined as being within 1KB of target gene, BUT expanding to 10kb didn't have much effect
eqtl.genes <- unique(eqtl.snp.table$GENE)
cis.eqtl.table <- eqtl.snp.table[row.names(eqtl.snp.table) %in% row.names(snp.gene.map),]
cis.snp.gene.map <- snp.gene.map[row.names(snp.gene.map) %in% row.names(cis.eqtl.table),]
cis.snp.gene.map <- cis.snp.gene.map[row.names(cis.eqtl.table),]
if(gene.mode) {
  if(!all(cis.eqtl.table$TARGET.DIST == cis.snp.gene.map$DIST)) {
    cat("WARNING: Some target distances computed from SNP and Gene coordinates do not match up with distances reported in", my.args$MAP, "\n")
    cat("Correcting to the distances reported in Map file.\n")
    cis.eqtl.table$TARGET.DIST <- cis.snp.gene.map$DIST
  }
} else {
  # Non-Gene Mode (shouldn't have any cis eQTLs unless there is a problem with map file)
  # Just pull TARGET.DIST from the map
  cis.eqtl.table$TARGET.DIST <- cis.snp.gene.map$DIST
}


trans.eqtl.table <- eqtl.snp.table[setdiff(row.names(eqtl.snp.table),row.names(cis.eqtl.table)),]
cat("\n")
cat(nrow(cis.eqtl.table), " of ", nrow(eqtl.snp.table), " eQTL SNPs (", round(nrow(cis.eqtl.table)*100/nrow(eqtl.snp.table)), "%) are cis based on map in", my.args$MAP, "\n", sep="")
cat(length(unique(cis.eqtl.table$GENE)), " of ", length(eqtl.genes), " gene features (", round(length(unique(cis.eqtl.table$GENE))*100/length(eqtl.genes)), "%) are linked to at least one cis eQTL.\n\n", sep="")
cat(nrow(trans.eqtl.table), " of ", nrow(eqtl.snp.table), " eQTL SNPs (", round(nrow(trans.eqtl.table)*100/nrow(eqtl.snp.table)), "%) are trans.\n", sep="")
cat(length(unique(trans.eqtl.table$GENE)), " of ", length(eqtl.genes), " gene features (", round(length(unique(trans.eqtl.table$GENE))*100/length(eqtl.genes)), "%) are linked to at least one trans eQTL.\n", sep="")
if(gene.mode) {
  cat(" ", sum(trans.eqtl.table$TARGET.DIST < Inf), " of ", nrow(trans.eqtl.table), " trans eQTL SNPs (", round(mean(trans.eqtl.table$TARGET.DIST < Inf)*100), "%) are on same chromosome as target gene.\n", sep="")
  cat(" Median distance to target gene =", median(trans.eqtl.table[trans.eqtl.table$TARGET.DIST < Inf,"TARGET.DIST"]), "bp\n\n")
}
cat(length(intersect(unique(cis.eqtl.table$GENE), trans.eqtl.table$GENE)), " of ", length(eqtl.genes), " genes (", round(length(intersect(unique(cis.eqtl.table$GENE), trans.eqtl.table$GENE))*100/length(eqtl.genes)), "%) have both cis and trans eQTL.\n", sep="")
cat(length(setdiff(unique(cis.eqtl.table$GENE), trans.eqtl.table$GENE)), " of ", length(eqtl.genes), " genes (", round(length(setdiff(unique(cis.eqtl.table$GENE), trans.eqtl.table$GENE))*100/length(eqtl.genes)), "%) have cis eQTL ONLY.\n", sep="")
cat(length(setdiff(unique(trans.eqtl.table$GENE), cis.eqtl.table$GENE)), " of ", length(eqtl.genes), " genes (", round(length(setdiff(unique(trans.eqtl.table$GENE), cis.eqtl.table$GENE))*100/length(eqtl.genes)), "%) have trans eQTL ONLY.\n\n", sep="")
rm(eqtl.snp.table)

# Split distal eQTLs (no linked gene) out of trans.eqtl.table
trans.eqtl.total <- nrow(trans.eqtl.table)
distal.eqtl.table <- trans.eqtl.table[!(trans.eqtl.table$SNP %in% snp.gene.map$SNP),]
trans.eqtl.table <- trans.eqtl.table[trans.eqtl.table$SNP %in% snp.gene.map$SNP,]

cat(nrow(trans.eqtl.table), " of ", trans.eqtl.total, " trans eQTLs (", round(nrow(trans.eqtl.table)*100/trans.eqtl.total), "%) are linked to at least one candidate gene.\n", sep="")
cat(nrow(distal.eqtl.table), " of ", trans.eqtl.total, " trans eQTLs (", round(nrow(distal.eqtl.table)*100/trans.eqtl.total), "%) are not associated with any genes.\n", sep="")
if(gene.mode) {
  cat(" ", sum(distal.eqtl.table$TARGET.DIST < Inf), " (", round(mean(distal.eqtl.table$TARGET.DIST < Inf)*100), "%) are on the same chromosome as target gene.\n", sep="")
  cat(" Median distance for these = ", median(distal.eqtl.table[distal.eqtl.table$TARGET.DIST < Inf,"TARGET.DIST"]), " bp\n", sep="")
  for(d in 4:6) {
    cat(" ", sum(distal.eqtl.table$TARGET.DIST < (10^d)), " (", round(mean(distal.eqtl.table$TARGET.DIST < (10^d))*100), "%) are within ", (10^(d-3)), "kb of target gene.\n", sep="")
  }
  cat("\n")
}


# Now map upstream genes for all remaining trans eQTLs

# Get the list of trans eQTL SNP IDs
trans.snp.table <- data.frame(row.names=unique(trans.eqtl.table$SNP))

# TEMP TESTING:
# trans.snp.table <- data.frame(row.names=unique(trans.eqtl.table$SNP)[1:100])
# trans.eqtl.table <- trans.eqtl.table[trans.eqtl.table$SNP %in% row.names(trans.snp.table),]

# DEPRECATED:
# snp.split <- strsplit(row.names(trans.snp.table), split="_", fixed=T)
# trans.snp.table$CHR <- unlist(lapply(snp.split, "[", 1))
# trans.snp.table$POS <- as.integer(unlist(lapply(snp.split, "[", 2)))
# 
# # determine all genes within 1kb for all SNPs
# cat("Mapping all trans eQTLs to all genes within", my.args$DIST, "bp...\n")
# snp.gene.map <- foreach(snp=row.names(trans.snp.table)) %dopar% {
#   snp.chr.genes <- gene.bodies[gene.bodies$CHR==trans.snp.table[snp,"CHR"],]
#   snp.pos <- trans.snp.table[snp,"POS"]
#   # Compute distance from SNP to all genes on same CHR
#   snp.chr.genes[snp.chr.genes$START > snp.pos,"DIST"] <- snp.chr.genes[snp.chr.genes$START > snp.pos,"START"] - snp.pos
#   snp.chr.genes[snp.chr.genes$END < snp.pos, "DIST"] <- snp.pos - snp.chr.genes[snp.chr.genes$END < snp.pos,"END"]
#   snp.chr.genes[(snp.chr.genes$START <= snp.pos) & (snp.chr.genes$END >= snp.pos),"DIST"] <- 0
#   stopifnot(sum(is.na(snp.chr.genes$DIST))==0)
#   snp.chr.genes <- snp.chr.genes[snp.chr.genes$DIST <= my.args$DIST,]
#   # Now further annotate with part of gene for each case
#   if(nrow(snp.chr.genes) > 0) {
#     snp.chr.genes$PART <- as.character(NA)
#     overlap.exons <- unlist(lapply(row.names(snp.chr.genes), function(g){
#       my.exons <- gene.exons[[g]]
#       any((my.exons$START <= snp.pos) & (my.exons$END >= snp.pos))
#     }))
#     snp.chr.genes[overlap.exons,"PART"] <- "EXON"
#     snp.chr.genes[(snp.chr.genes$DIST==0) & !overlap.exons,"PART"] <- "INTRON"
#     distal.genes <- snp.chr.genes$DIST > 0
#     upstream.genes <- distal.genes & (((snp.chr.genes$STRAND=="+") & (snp.chr.genes$START > snp.pos)) | ((snp.chr.genes$STRAND=="-") & (snp.chr.genes$END < snp.pos)))
#     snp.chr.genes[upstream.genes,"PART"] <- "5PRIME"
#     snp.chr.genes[is.na(snp.chr.genes$PART),"PART"] <- "3PRIME"
#   } else {
#     snp.chr.genes[,"PART"] <- character()
#   }
#   return(snp.chr.genes)
# }
# names(snp.gene.map) <- row.names(trans.snp.table)
# 
# # TO DO: Fold this into loop above?
# cat("Identifying proximal genes where the trans eQTL is also a cis eQTL...\n")
# snp.gene.map <- foreach(snp=names(snp.gene.map)) %dopar% {
#   snp.genes <- snp.gene.map[[snp]]
#   # Is this SNP a cis eQTL for any of these genes?
#   snp.cis.eqtls <- cis.eqtl.table[cis.eqtl.table$SNP==snp,]
#   snp.genes[,"CIS.EQTL"] <- row.names(snp.genes) %in% snp.cis.eqtls$GENE
#   return(snp.genes)
# }
# names(snp.gene.map) <- row.names(trans.snp.table)


# NEW: Use snp.gene.map to map trans eQTLs to nearest genes
cat("Mapping", nrow(trans.snp.table), "trans eQTLs to candidate genes based on",my.args$MAP,"and cis eQTLs...\n")
start.time <- proc.time()
trans.eqtl.gene.map <- foreach(snp=row.names(trans.snp.table)) %dopar% {
  snp.genes <- snp.gene.map[snp.gene.map$SNP==snp,2:5]
  row.names(snp.genes) <- snp.genes$GENE
  snp.cis.eqtls <- cis.eqtl.table[cis.eqtl.table$SNP==snp,]
  snp.genes[,"CIS.EQTL"] <- row.names(snp.genes) %in% snp.cis.eqtls$GENE
  return(snp.genes)
}
names(trans.eqtl.gene.map) <- row.names(trans.snp.table)

# TO DO: Fold into loops above?
# Now go through and pick a best gene (or genes) for each SNP
# Simplifying the rules heavily here - if SNP is eQTL for any of the underlying genes, use only those, otherwise use everything else!
trans.snp.genes <- foreach(snp=row.names(trans.snp.table), .combine='rbind') %dopar% {
  snp.genes <- trans.eqtl.gene.map[[snp]]
  snp.best.gene <- data.frame(SNP.GENE="", SNP.GENE.REASON="", SNP.GENE.DIST="", SNP.GENE.PART="", SNP.GENE.STRAND="", row.names=snp)
  if(nrow(snp.genes) > 0) {
    if(any(snp.genes$CIS.EQTL)) {
      snp.genes <- snp.genes[snp.genes$CIS.EQTL,]
      snp.gene.reason <- "CIS"
    } else {
      snp.gene.reason <- "PROX"
    }
    snp.best.gene[,"SNP.GENE"] <- paste(row.names(snp.genes), collapse=",")
    snp.best.gene[,"SNP.GENE.REASON"] <- paste(rep(snp.gene.reason, times=nrow(snp.genes)), collapse=",")
    snp.best.gene[,"SNP.GENE.DIST"] <- paste(snp.genes$DIST, collapse=",")
    snp.best.gene[,"SNP.GENE.PART"] <- paste(snp.genes$GENE.PART, collapse=",")
    snp.best.gene[,"SNP.GENE.STRAND"] <- paste(snp.genes$GENE.STRAND, collapse=",")
  }
  return(snp.best.gene)
}
end.time <- proc.time()
cat("Finished in:\n")
print(end.time - start.time)
cat("\n")

# Now map these over to trans.eqtl.table - Should ONLY be mappable SNPs now
trans.eqtl.table <- cbind(trans.eqtl.table, trans.snp.genes[trans.eqtl.table$SNP,])
stopifnot(sum(is.na(trans.eqtl.table$SNP.GENE))==0)
stopifnot(sum(trans.eqtl.table$SNP.GENE=="")==0)

# DEPRECATED
# # Split into Gene-Gene connections ("linked") from everything else ("orphan")
# linked.trans.eqtl.table <- trans.eqtl.table[trans.eqtl.table$SNP.GENE != "",]
# orphan.trans.eqtl.table <- trans.eqtl.table[trans.eqtl.table$SNP.GENE == "",1:10]

cat(sum(!grepl(",", trans.eqtl.table$SNP.GENE)), " of ", nrow(trans.eqtl.table), " trans eQTLs (", round(mean(!grepl(",", trans.eqtl.table$SNP.GENE))*100), "%) are linked to exactly one gene.\n", sep="")
cat(sum(grepl(",", trans.eqtl.table$SNP.GENE)), " of ", nrow(trans.eqtl.table), " trans eQTLs (", round(mean(grepl(",", trans.eqtl.table$SNP.GENE))*100), "%) are linked to multiple candidate genes.\n", sep="")
cat(sum(grepl("CIS", trans.eqtl.table$SNP.GENE.REASON)), " of ", nrow(trans.eqtl.table), " trans eQTLs (", round(mean(grepl("CIS", trans.eqtl.table$SNP.GENE.REASON))*100), "%) are linked based on cis eQTL target.\n", sep="")
cat(sum(!grepl("CIS", trans.eqtl.table$SNP.GENE.REASON)), " of ", nrow(trans.eqtl.table), " trans eQTLs (", round(mean(!grepl("CIS", trans.eqtl.table$SNP.GENE.REASON))*100), "%) are linked based on proximity.\n\n", sep="")

# I'm quite surprised at how few cis eQTLs I'm finding, maybe the first 100 genes are just not representative, I don't know yet
# I suspect there are at least some cases where the trans eQTL will be linked to a nearby novel gene
# that is actually just a novel start exon for same gene that wasn't properly merged in transcriptome build
#  - should try to estimate how many cases look like this...
# TO DO: Dump some of the more useful cis eQTL columns to a table (NEED OUTPUT PARAMS)

# NOTE: PART OF THE REASON FOR THIS MIGHT BE THE GENES ON NON-STD CHROMOSOMES!!
# MAYBE WE SHOULD DROP THEM FROM ANALYSIS OR AT LEAST COMPUTE % OF GENES ON STD CHROMS FOR % eQTL

# TO DO: It would be useful to compute where in/around target gene the cis eQTLs fall?

# Create output table for cis eQTLs
cat("Contructing final output table for cis eQTLs...\n")
cis.eqtl.output <- data.frame(GENE=unique(cis.eqtl.table$GENE))
if(gene.mode) {
  cis.eqtl.output <- cbind(cis.eqtl.output, gene.bodies[cis.eqtl.output$GENE,])
}
cis.eqtl.output[,"CIS.SNPS"] <- foreach(gene=cis.eqtl.output$GENE, .combine='c') %dopar% {
  paste(cis.eqtl.table[cis.eqtl.table$GENE==gene,"SNP"], collapse=",")
}
cis.eqtl.output[,"SNP.PVALS"] <- foreach(gene=cis.eqtl.output$GENE, .combine='c') %dopar% {
  paste(cis.eqtl.table[cis.eqtl.table$GENE==gene,"PVAL"], collapse=",")
}
if(gene.mode) {
  cis.eqtl.output[,"SNP.DISTS"] <- foreach(gene=cis.eqtl.output$GENE, .combine='c') %dopar% {
    paste(cis.eqtl.table[cis.eqtl.table$GENE==gene,"TARGET.DIST"], collapse=",")
  }
}
cis.eqtl.file <- paste0(my.args$OUTPUT, ".cis.eqtls.txt")
cat("Writing cis eQTL results to:", cis.eqtl.file, "\n")
write.table(cis.eqtl.output, cis.eqtl.file, sep="\t", row.names=F, quote=F)
cat("\n")

# Output table trans eQTLs (linked to other genes) (leave in current format?)
trans.eqtl.file <- paste0(my.args$OUTPUT, ".trans.eqtls.txt")
cat("Writing trans eQTL results to:", trans.eqtl.file, "\n")
write.table(trans.eqtl.table, trans.eqtl.file, sep="\t", row.names=F, quote=F)
cat("\n")

# TO DO: Create simpler table of the gene-gene network inferred from eQTLs

# TO DO: Output table of distal (orphan.trans) eQTLs (pack in similar format to cis eQTLs)
cat("Contructing final output table for distal eQTLs...\n")
distal.eqtl.output <- data.frame(GENE=unique(distal.eqtl.table$GENE))
if(gene.mode) {
  distal.eqtl.output <- cbind(distal.eqtl.output, gene.bodies[distal.eqtl.output$GENE,])
}
distal.eqtl.output[,"DISTAL.SNPS"] <- foreach(gene=distal.eqtl.output$GENE, .combine='c') %dopar% {
  paste(distal.eqtl.table[distal.eqtl.table$GENE==gene,"SNP"], collapse=",")
}
distal.eqtl.output[,"SNP.PVALS"] <- foreach(gene=distal.eqtl.output$GENE, .combine='c') %dopar% {
  paste(distal.eqtl.table[distal.eqtl.table$GENE==gene,"PVAL"], collapse=",")
}
if(gene.mode) {
  distal.eqtl.output[,"SNP.DISTS"] <- foreach(gene=distal.eqtl.output$GENE, .combine='c') %dopar% {
    paste(distal.eqtl.table[distal.eqtl.table$GENE==gene,"TARGET.DIST"], collapse=",")
  }
}
distal.eqtl.file <- paste0(my.args$OUTPUT, ".distal.eqtls.txt")
cat("Writing distal eQTL results to:", distal.eqtl.file, "\n")
write.table(distal.eqtl.output, distal.eqtl.file, sep="\t", row.names=F, quote=F)
cat("\n")

cat("\n\nScript completed successfully!\n\n")
print(proc.time())
