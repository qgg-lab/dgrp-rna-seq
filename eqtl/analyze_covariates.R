#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 2/21/17
#
# plot_mmc_results.R
#
# Goal: Load line means, clustering results, and analyze covariates of clustering and other aspects
# 
# Usage:
# Rscript plot_mmc_results.R EXPR=subdir/data_line_means.txt [OPTIONS]
#  EXPR=      Path to file of line means
#  H2FILE=    Path to H2 model results file
#  CLUST=     Path to MMC/WGCNA output file in CSV format (has cluster assignments)
#  PCA=       Path to PCA results (the gene results containing Perc.Var)
#  SEX=       Which sex in H2 to use (inferred from EXPR file name by default)
#  GENES=   path to file of gene info, must contain LENGTH column for FPKM normalization (if row IDs don't match up completely, skip FPKM normalization)
#  TISSUES= Path to Fly Atlas expression data, defaults to Tissue_Deconvolution/FlyAtlas/FlyAtlas_gene_bodyrle_table.txt
#  

# TO DO: Plot various properties across clusters and vs cluster size, and vs cluster properties

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

options(stringsAsFactors=F)

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript analyze_covariates.R EXPR=subdir/data_line_means.txt CLUST=subdir/mmc/data_mmc.csv"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_F_line_means.txt","H2FILE=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_model_results.txt","CLUST=known_all_novel_genes/mmc/combined_samples_known_novel_fpkm_VR_WolAdj_F_mmc.csv","PCA=known_all_novel_genes/pca/combined_samples_known_novel_fpkm_VR_WolAdj_F_PC_gene_results.txt","GENES=~/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-transcriptome-r5.57-gene-info.txt")
# my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_F_line_means.txt","H2FILE=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_model_results.txt","CLUST=known_all_novel_genes/wgcna/combined_samples_known_novel_fpkm_VR_WolAdj_F_wgcna.csv","GENES=~/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-transcriptome-r5.57-gene-info.txt")
if(length(my.args) < 2) {
  stop("Requires at least two arguments. ", usageStr)
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

# Make sure EXPR param exists
if(!("EXPR" %in% names(my.args))) {
  stop("Missing EXPR parameter. ", usageStr)
}

#  H2FILE=  The .._model_results.txt file to use for H2 and FDR thresholds
setDefault("H2FILE","")

#  SEX=   M/F/Pooled - tells which columns in H2FILE to use.  Tries to infer from line_means.txt file if not specified.
inferredSex <- regmatches(my.args$EXPR, regexec(".*[_]([A-Za-z]+)[_]line[_]means[.]txt", my.args$EXPR))[[1]][2]
if(inferredSex == "SexAvg") {inferredSex <- "Pooled"}
if(!(inferredSex %in% c("M","F","Pooled"))) {inferredSex <- ""}
# Change SexAvg to Pooled, set to "" if not one of the main values
setDefault("SEX",inferredSex)

# Default CLUST file is empty (no cluster analysis)
setDefault("CLUST", "")

setDefault("PCA","")

setDefault("GENES", "~ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-transcriptome-r5.57-gene-info.txt")

setDefault("TISSUES", "Tissue_Deconvolution/FlyAtlas/FlyAtlas_gene_bodyrle_table.txt")


# --- LOAD INPUT FILES --- #

# Load the line means
cat("Loading line means from:", my.args$EXPR, "\n")
expr.table <- read.table(my.args$EXPR, header=T, sep="\t", row.names=1, as.is=T)
colnames(expr.table) <- sub("^X","", colnames(expr.table))
cat("Loaded line means for", nrow(expr.table), "genes across", ncol(expr.table), "samples.\n\n")

if((my.args$H2FILE != "") & file.exists(my.args$H2FILE)) {
  cat("Loading model results from:", my.args$H2FILE, "\n")
  h2.table <- read.table(my.args$H2FILE, header=T, sep="\t", row.names=1)
  h2.sex.cols <- grep(paste0("^",my.args$SEX,"[.]"), colnames(h2.table), value=T)
  cat("Limiting to", length(h2.sex.cols), "columns from", my.args$SEX, "model.\n")
  h2.table <- h2.table[,h2.sex.cols]
  colnames(h2.table) <- sub(paste0("^",my.args$SEX,"[.]"), "", colnames(h2.table))
} else {
  cat("WARNING: H2FILE parameter not specified or file not found, skipping comparison to gen var model results.\n")
  h2.table <- NULL
}

if((my.args$GENES != "") & file.exists(my.args$GENES)) {
  feature.table <- read.table(my.args$GENES, header=T, sep="\t", row.names=1)
  cat("Loaded feature length for",nrow(feature.table),"genes/features from",my.args$GENES,"\n")
} else {
  cat("WARNING: GENES parameter not specified or file not found.\n")
  feature.table <- NULL
}

if((my.args$TISSUES != "") & file.exists(my.args$TISSUES)) {
  tissue.table <- read.table(my.args$TISSUES, header=T, sep="\t", row.names=1)
  # Drop everything except OveryMean, TestisMean, AccMean
  tissue.table <- tissue.table[,c("OvaryMean","TestisMean","AccMean")]
  # Split any row.name that represents multiple genes
  tissue.table.multi <- tissue.table[grepl(" \\/\\/\\/ ", row.names(tissue.table)),]
  tissue.table <- tissue.table[!grepl(" \\/\\/\\/ ", row.names(tissue.table)),]
  gene.multi.map <- strsplit(row.names(tissue.table.multi), split=" /// ", fixed=T)
  names(gene.multi.map) <- row.names(tissue.table.multi)
  tissue.multi.rowcount <- unlist(lapply(gene.multi.map, length))
  tissue.table.split <- tissue.table.multi[rep(names(tissue.multi.rowcount), times=tissue.multi.rowcount),]
  tissue.table.split[,"FBGN"] <- unlist(gene.multi.map)
  # Remove all genes that are redundant/ambiguous
  tissue.table.split <- tissue.table.split[!duplicated(tissue.table.split$FBGN),]
  tissue.table.split <- tissue.table.split[!(tissue.table.split$FBGN %in% row.names(tissue.table)),]
  row.names(tissue.table.split) <- tissue.table.split$FBGN
  tissue.table <- rbind(tissue.table, tissue.table.split[,colnames(tissue.table)])
  cat("Loaded tissue-specificity data for", nrow(tissue.table), "genes from", my.args$TISSUES, "\n")
} else {
  cat("WARNING: TISSUES parameter not specified or file not found, skipping comparison to tissue-specificity.\n")
  tissue.table <- NULL
}

if((my.args$CLUST != "") & file.exists(my.args$CLUST)) {
  clust.table <- read.table(my.args$CLUST, header=T, sep=",", row.names=1)
  # Compute cluster size
  clust.table[,"Cluster.Size"] <- table(clust.table$Module)[as.character(clust.table$Module)]
  clust.table[clust.table$Module==0,"Cluster.Size"] <- 0
  cat("Loaded cluster assignments for", nrow(clust.table), "genes from", my.args$CLUST, "\n")
} else {
  cat("WARNING: CLUST parameter not specified or file not found, skipping cluster-wise analysis.\n")
  clust.table <- NULL
}

if((my.args$PCA != "") & file.exists(my.args$PCA)) {
  pca.table <- read.table(my.args$PCA, header=T, sep="\t", row.names=1)
  cat("Loaded PCA results for", nrow(pca.table), "genes from", my.args$PCA, "\n")
} else {
  cat("WARNING: PCA parameter not specified or file not found, skipping PC1 analysis.\n")
  pca.table <- NULL
}

cat("\n")


# --- ANALYZE COVARIATES --- #

# Put flags in seperate data structure
gene.flags <- expr.table$FLAG
names(gene.flags) <- row.names(expr.table)
expr.table <- expr.table[,setdiff(colnames(expr.table),"FLAG")]

# Compute mean expr and variance for each gene
gene.means <- apply(expr.table, 1, mean)
gene.var <- apply(expr.table, 1, var)

### DEPRECATED ###
# # Compute the variance
# cat("Performing PCA analysis...\n")
# expr.prcomp <- prcomp(t(expr.table), scale=T, retx=T)
# expr.pc.var <- expr.prcomp$sdev ** 2
# expr.pc.pvar <- expr.pc.var / sum(expr.pc.var)
# expr.pc.cvar <- unlist(lapply(1:length(expr.pc.pvar), function(i){sum(expr.pc.pvar[1:i])}))
# # expr.pc.table <- t(expr.prcomp$x)
# 
# expr.pc.output <- cbind(data.frame(PC=colnames(expr.prcomp$x), Var=expr.pc.var, PropVar=expr.pc.pvar, CumPropVar=expr.pc.cvar))
# cat("Top PCs:\n")
# print(head(expr.pc.output[,1:4]))
# 
# # Extract first PC
# expr.pc1 <- expr.prcomp$x[,1]
# 
# # TO DO: This part is still pretty slow, can either parallelize or try to find other ways to speed it up
# cat("Regressing PC1 out of", nrow(expr.table), "profiles")
# expr.pc1.lm <- foreach(gene=row.names(expr.table), .combine=rbind) %do% {
#   pc1.lm <- lm(t(expr.table)[,gene] ~ expr.pc1)
#   pc1.perc.var <- (gene.var[gene] - var(pc1.lm$residuals)) / gene.var[gene]
#   pc1.pval <- anova(pc1.lm)[["Pr(>F)"]][1]
#   c(pc1.perc.var, pc1.pval)
# }
# colnames(expr.pc1.lm) <- c("PC1.Perc.Var", "PC1.PVal")
# row.names(expr.pc1.lm) <- row.names(expr.table)
# expr.pc1.lm <- as.data.frame(expr.pc1.lm)
# 
# cat(sum(expr.pc1.lm$PC1.PVal <= 0.01), "genes have significant PC1 component (P < 0.01)\n")
# cat("Distribution of % Variance Explained:\n")
# cat(min(expr.pc1.lm$PC1.Perc.Var), boxplot.stats(expr.pc1.lm$PC1.Perc.Var)$stats, max(expr.pc1.lm$PC1.Perc.Var), "\n")

# Make a table of features to correlate against each other:
gene.metrics <- data.frame(EXPR.MEAN=gene.means, EXPR.VAR=gene.var)

if(!is.null(pca.table)) {
  gene.metrics <- cbind(pca.table[row.names(gene.metrics),c("PC1.Perc.Var","PC2.Perc.Var")], gene.metrics)
}

if(!is.null(feature.table)) {
  gene.metrics <- cbind(gene.metrics, LENGTH=feature.table[row.names(gene.metrics),"LENGTH"])
}

if(!is.null(h2.table)) {
  perc.var.cols <- grep("[.]Perc[.]Var$", colnames(h2.table), value=T)
  if("H2" %in% colnames(h2.table)) {
    perc.var.cols <- c(perc.var.cols, "H2")
  }
  gene.metrics <- cbind(gene.metrics, h2.table[row.names(gene.metrics),perc.var.cols])
}

if(!is.null(tissue.table)) {
  gene.metrics <- cbind(gene.metrics, tissue.table[row.names(gene.metrics),])
}

if(!is.null(clust.table)) {
  gene.metrics <- cbind(gene.metrics, clust.table[row.names(gene.metrics),c("Average.Degree","Degree","Cluster.Size")])
}

cat("Pearson correlations between gene properties:\n")
print(cor(gene.metrics, use="pairwise.complete.obs"))
cat("\n\n")

cat("Spearman correlations between gene properties:\n")
cor(gene.metrics, method="spearman", use="complete.obs")
cat("\n\n")

# Put together table of -log10 PVals
# TO DO: FlyAtlas P-values would be useful here, but more complicated to pull in
gene.signifs <- data.frame(row.names=row.names(expr.table))

if(!is.null(pca.table)) {
  gene.signifs[,"PC1.PVal"] <- pca.table[row.names(expr.table),"PC1.PVal"]
  gene.signifs[,"PC2.PVal"] <- pca.table[row.names(expr.table),"PC2.PVal"]
}

if(!is.null(h2.table)) {
  pval.cols <- grep("[.]PVal$", colnames(h2.table), value=T)
  gene.signifs <- cbind(gene.signifs, h2.table[row.names(gene.signifs),pval.cols])
}

for(j in 1:ncol(gene.signifs)) {
  if(any(gene.signifs[,j] == 0)) {
    j.zero <- gene.signifs[,j] == 0
    j.min.nonzero <- min(gene.signifs[!j.zero,j])
    gene.signifs[j.zero,j] <- j.min.nonzero/10
  }
}

gene.signifs <- -1 * log10(gene.signifs)

if(ncol(gene.signifs) > 1) {
  cat("Pearson correlations (-log10 P-values):\n")
  print(cor(gene.signifs))
  cat("\n\n")
  
  cat("Spearman correlations (-log10 P-values):\n")
  print(cor(gene.signifs, method="spearman"))
  cat("\n\n")
}


# Put together table of categories
gene.categories <- data.frame(FLAG=gene.flags, row.names=names(gene.flags))

if(!is.null(feature.table)) {
  gene.categories[,"CLASS"] <- feature.table[row.names(gene.categories),"CLASS"]
}

if(!is.null(tissue.table)) {
  # Everything >2x enriched is tissue specific for now
  tissue.groups <- tissue.table >= 1
  colnames(tissue.groups) <- sub("Mean$", "Specific", colnames(tissue.groups))
  for(j in colnames(tissue.groups)) {
    j.yes <- intersect(row.names(tissue.groups)[tissue.groups[,j]], row.names(gene.categories))
    j.no <- intersect(row.names(tissue.groups)[!tissue.groups[,j]], row.names(gene.categories))
    gene.categories[,j] <- "ABSENT"
    gene.categories[j.no,j] <- "NO"
    gene.categories[j.yes,j] <- "YES"
  }
}

if(!is.null(clust.table)) {
  gene.categories[,"CLUSTER"] <- "FILTERED"
  for(k in sort(unique(clust.table$Module))) {
    k.genes <- intersect(row.names(clust.table)[clust.table$Module==k], row.names(gene.categories))
    gene.categories[k.genes,"CLUSTER"] <- paste0("CLUST", k)
  }
}

# For each category, look at distribution of each metric and significance value
# TO DO: Could plot this in one big ream of PDFs?
cat("\n\n")
cat("--------------------------------------------------------------------------\n")
cat("Comparing all gene categorizations vs all metrics and significance values.\n")
cat("--------------------------------------------------------------------------\n\n")

for(j in colnames(gene.categories)) {
  if(length(unique(gene.categories[,j]))==1) {
    cat("SKIPPING", j, " - does NOT contain multiple category values.\n\n")
  } else {
    cat("Comparing metrics by",j,"column with groupings:\n")
    j.groups <- sort(unique(gene.categories[,j]))
    j.genes <- lapply(j.groups, function(grp){
      row.names(gene.categories)[gene.categories[,j]==grp]
    })
    names(j.genes) <- j.groups
    print(unlist(lapply(j.genes, length)))
    cat("\n")
    
    for(metric in colnames(gene.metrics)) {
      j.metric.values <- lapply(j.genes, function(x){gene.metrics[x,metric]})
      j.metric.valid <- unlist(lapply(j.metric.values, function(x){sum(!is.na(x))}))
      # Can only apply to cases where there are at least 2 valid metric values in at least 2 groups
      if(sum(j.metric.valid > 1) > 1) {
        j.metric.pval <- as.data.frame(anova(lm(gene.metrics[,metric] ~ gene.categories[,j])))[1,5]
        if(j.metric.pval <= 0.001) {
          cat("Significant association between", j, "and", metric, "(P =", j.metric.pval, ")\n")
          
          j.metric.means <- unlist(lapply(j.metric.values, mean, na.rm=T))
          j.metric.medians <- unlist(lapply(j.metric.values, median, na.rm=T))
          j.metric.vars <- unlist(lapply(j.metric.values, var, na.rm=T))
          cat("Group means:\n")
          print(j.metric.means)
          cat("Group medians:\n")
          print(j.metric.medians)
          cat("Group variance:\n")
          print(j.metric.vars)
          cat("\n")
        } else {
          cat("No significant associate between", j, "and", metric, "(P =", j.metric.pval, ")\n\n")
        }
      } else {
        cat("Insufficient values to test", j , "vs", metric, "\n\n")
      }
    }
  }
  cat("\n\n")
}


cat("\n\nScript completed successfully!\n\n")
print(proc.time())
cat("\n")
