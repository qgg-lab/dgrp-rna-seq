#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 6/28/16
#
# plot_mmc_results.R
#
# Goal: Load MMC clustering results, draw heatmap of clusters
# 
# Usage:
# Rscript plot_mmc_results.R EXPR=subdir/data_line_means.txt MMC=subdir/mmc/data_mmc.csv
#  EXPR=    Path to file of expression data used as input for MMC
#  MMC=     Path to MMC output file in CSV format (has cluster assignments)
#  

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

# library(RColorBrewer)
library(colorRamps)

options(stringsAsFactors=F)

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript plot_mmc_results.R EXPR=subdir/data_line_means.txt MMC=subdir/mmc/data_mmc.csv"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_F_line_means.txt","MMC=known_all_novel_genes/mmc/combined_samples_known_novel_fpkm_VR_WolAdj_F_mmc.csv")
# my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_M_line_means.txt","MMC=known_all_novel_genes/wgcna/combined_samples_known_novel_fpkm_VR_WolAdj_M_wgcna.csv")
# my.args <- c("EXPR=microbiome_species_filtered/genVar/combined_microbe_filtered_rpm_WolAdj_SexAvg_line_means.txt","MMC=microbiome_species_filtered/mmc/combined_microbe_filtered_rpm_WolAdj_SexAvg_mmc.csv")
# my.args <- c("EXPR=microbiome_species_filtered/genVar/combined_microbe_filtered_rpm_WolAdj_SexAvg_line_means.txt","MMC=microbiome_species_filtered/wgcna/combined_microbe_filtered_rpm_WolAdj_SexAvg_wgcna.csv")

if(length(my.args) < 2) {
  stop("Requires at least two arguments. ", usageStr)
}

# Parse param names and values
argSplit=strsplit(my.args, split="=", fixed=T)
my.args <- lapply(argSplit, "[", 2)
names(my.args) <- unlist(lapply(argSplit, "[", 1))

# Make sure EXPR param exists
if(!("EXPR" %in% names(my.args))) {
  stop("Missing EXPR parameter. ", usageStr)
}

if(!("MMC" %in% names(my.args))) {
  stop("Missing MMC parameter. ", usageStr)
}

# Output files based on input names
my.args$SORTFILE=sub("[.]csv$", "_sorted_heatmap.png", my.args$MMC)
my.args$SMOOTHFILE=sub("[.]csv$", "_smoothed_heatmap.png", my.args$MMC)
my.args$CROPSORTFILE=sub("[.]csv$", "_cropped_sorted_heatmap.png", my.args$MMC)
my.args$CROPSMOOTHFILE=sub("[.]csv$", "_cropped_smoothed_heatmap.png", my.args$MMC)
my.args$PALFILE=sub("[.]csv$", "_palette.png", my.args$MMC)
my.args$CORFILE=sub("[.]csv$", "_within_module_cor.pdf", my.args$MMC)

# --- LOAD INPUT FILES --- #

# Load the line means
cat("Loading line means from:", my.args$EXPR, "\n")
expr.table <- read.table(my.args$EXPR, header=T, sep="\t", row.names=1)
cat("Loaded line means for", nrow(expr.table), "genes across", ncol(expr.table), "samples.\n\n")
# Do the same transformations on row names to make sure things match up
# (This really only matters for microbiome)
row.names(expr.table) <- gsub(" ", ".", row.names(expr.table))
row.names(expr.table) <- gsub(",", "|", row.names(expr.table))

# Load the MMC cluster results
cat("Loading MMC cluster assignments:", my.args$MMC, "\n")
mmc.table <- read.csv(my.args$MMC, row.names=1, header=T)
cat("Loaded MMC cluster assignments for", nrow(mmc.table), "genes.\n\n")
row.names(mmc.table) <- gsub(" ", ".", row.names(mmc.table))
row.names(mmc.table) <- gsub(",", "|", row.names(mmc.table))

# Drop the genes missing from either table - keep order the same as in MMC table
common.genes <- intersect(row.names(mmc.table), row.names(expr.table))
expr.drop <- length(setdiff(row.names(expr.table), common.genes))
mmc.drop <- length(setdiff(row.names(mmc.table), common.genes))

if(expr.drop > 0) {
  cat(expr.drop, "genes have line means but are not in MMC results, they will be dropped.\n")
}
if(mmc.drop > 0) {
  cat("WARNING:", mmc.drop, "genes are included in MMC results but do not have line means, are you sure the MMC results came from these line means?\n")
}
cat("\n")

expr.table <- expr.table[common.genes,]
mmc.table <- mmc.table[common.genes,]

# If expr.table contains a FLAG column, report the flag counts then drop it
if("FLAG" %in% colnames(expr.table)) {
  cat("Flag counts for clustered features:\n")
  print(table(expr.table$FLAG))
  cat("\n")
  expr.table <- expr.table[,setdiff(colnames(expr.table),"FLAG")]
}

# Confirm that we haven't disrupted intended MMC sorting
if(all(mmc.table$Module[mmc.table$Module > 0] == sort(mmc.table$Module[mmc.table$Module > 0], decreasing=F))) {
  cat("CHECK: Genes are already sorted by module.\n")
} else {
  cat("WARNING: Genes are not sorted by module!\n")
}

if(all(mmc.table$Average.Degree == sort(mmc.table$Average.Degree, decreasing=T))) {
  cat("CHECK: Clusters are sorted by avg degree.\n")
} else {
  cat("WARNING: Clusters are not sorted by avg degree!\n")
}

cluster.sorted <- T
for(k in unique(mmc.table$Module)) {
  cluster.degrees <- mmc.table[mmc.table$Module == k,"Degree"]
  cluster.sorted <- cluster.sorted & all(cluster.degrees == sort(cluster.degrees, decreasing=T))
}
if(cluster.sorted) {
  cat("CHECK: Genes are sorted by degree within clusters.\n")
} else {
  cat("WARNING: Genes are not sorted by degree within clusters!\n")
}

# Compute correlation matrix of all line means (currently using Pearson but this should be a command-line param?)
cat("Computing correlation table and reordering features within clusters...\n")
expr.cor <- cor(t(expr.table))

# k <- max(mmc.table$Module)

# TO DO: Try reordering modules by mean absolute correlation of all gene pairs?

# Loop over each module, compute mean of each gene against all others in module
# Reorder genes within module from strongest positive to strongest negative correlation
new.row.order <- c()
for(i in unique(mmc.table$Module)) {
  i.genes <- row.names(mmc.table)[mmc.table$Module==i]
  names(i.genes) <- i.genes
  i.avg.cor <- unlist(lapply(i.genes, function(x){
    i.other.genes <- setdiff(i.genes, x)
    mean(expr.cor[x,i.other.genes])
  }))
  new.row.order <- c(new.row.order, names(sort(i.avg.cor, decreasing=T)))
}
expr.cor <- expr.cor[new.row.order,new.row.order]
mmc.table <- mmc.table[new.row.order,]
# NOTE: Not reordering the expression table here, it doesn't get used again

# Draw a heatmap
# TO DO: Make a parameter for changing the palette choice
# heatmap.pal <- rev(colorRampPalette(brewer.pal(n=11, name="RdYlBu"))(100))
heatmap.pal <- matlab.like2(200)

# Old-style rainbow palette for module colors
# module.pal <- rainbow(n=k)
# module.pal <- module.pal[c(which((1:k %% 3)==1), which((1:k %% 3)==2), which((1:k %% 3)==0))]

# Instead just alternate grey and black
module.pal <- rep(c("lightgrey","black"), times=ceiling(length(unique(mmc.table$Module))/2))[1:length(unique(mmc.table$Module))]
names(module.pal) <- unique(mmc.table$Module)

cat("Drawing sorted heatmap to:", my.args$SORTFILE, "\n")
bitmap(my.args$SORTFILE, width=12, height=12, res=300)
heatmap(expr.cor, Rowv = NA, Colv=NA, scale="none", revC=T, zlim=c(-1,1), col=heatmap.pal, ColSideColors=module.pal[as.character(mmc.table$Module)], labRow=NA, labCol=NA)
dev.off()

# Smooth between module correlations
expr.cor.smooth <- expr.cor
modules <- unique(mmc.table$Module)
for(i in 1:(length(modules)-1)) {
  for(j in (i+1):length(modules)) {
    i.genes <- row.names(mmc.table)[mmc.table$Module==modules[i]]
    j.genes <- row.names(mmc.table)[mmc.table$Module==modules[j]]
    ij.avg <- mean(unlist(expr.cor[i.genes,j.genes]))
    expr.cor.smooth[i.genes,j.genes] <- ij.avg
    expr.cor.smooth[j.genes,i.genes] <- ij.avg
  }
}

cat("Drawing smoothed heatmap to:", my.args$SMOOTHFILE, "\n")
bitmap(my.args$SMOOTHFILE, width=12, height=12, res=300)
heatmap(expr.cor.smooth, Rowv = NA, Colv=NA, scale="none", revC=T, zlim=c(-1,1), col=heatmap.pal, ColSideColors=module.pal[as.character(mmc.table$Module)], labRow=NA, labCol=NA)
dev.off()

# Now crop to just clusters < 1000 genes and drop unsorted (cluster 0)
module.sizes <- table(mmc.table$Module)
crop.modules <- as.integer(names(module.sizes)[module.sizes < 1000])
crop.modules <- setdiff(crop.modules, 0)
cat("Cropping out the following modules:\n")
print(module.sizes[as.character(setdiff(mmc.table$Module,crop.modules))])
cat("\n")
crop.genes <- row.names(mmc.table)[mmc.table$Module %in% crop.modules]
crop.mmc.table <- mmc.table[crop.genes,]
# Renumber the modules:
crop.mod.map <- 1:length(unique(crop.mmc.table$Module))
names(crop.mod.map) <- unique(crop.mmc.table$Module)
crop.mmc.table$Module <- crop.mod.map[as.character(crop.mmc.table$Module)]
crop.expr.cor <- expr.cor[crop.genes,crop.genes]
crop.expr.cor.smooth <- expr.cor.smooth[crop.genes,crop.genes]

cat("Drawing cropped, sorted heatmap to:", my.args$CROPSORTFILE, "\n")
bitmap(my.args$CROPSORTFILE, width=12, height=12, res=300)
heatmap(crop.expr.cor, Rowv = NA, Colv=NA, scale="none", revC=T, zlim=c(-1,1), col=heatmap.pal, ColSideColors=module.pal[as.character(crop.mmc.table$Module)], labRow=NA, labCol=NA)
dev.off()

cat("Drawing cropped, smoothed heatmap to:", my.args$CROPSMOOTHFILE, "\n")
bitmap(my.args$CROPSMOOTHFILE, width=12, height=12, res=300)
heatmap(crop.expr.cor.smooth, Rowv = NA, Colv=NA, scale="none", revC=T, zlim=c(-1,1), col=heatmap.pal, ColSideColors=module.pal[as.character(crop.mmc.table$Module)], labRow=NA, labCol=NA)
dev.off()

# Copied this function from figPlot.R so that this script can be run by other users:
# drawHeatmapScale
# Function to draw the scale for a heatmap,
# given the zlim and color palette
drawHeatmapScale <- function(
  zlim,	# The extreme ends of the scale
  palette,	# The colors passed to the image or heatmap function
  xlab="Heatmap Scale",
  filename="",	# Set to a file name to draw to PNG
  height=0.6,
  width=2,
  cex=1.4,
  axis.limits=zlim,	# The extreme values for the axis (first two params of xaxp)
  axis.intervals=2	# Used as the third parameter for xaxp, number of intervals on the axis
) {
  # First, create a series of values with length of palette spanning zlim
  stepSz <- (zlim[2] - zlim[1]) / length(palette)
  
  gridlines <- zlim[1]
  for(i in 1:length(palette)) {
    gridlines[i+1] <- gridlines[i] + stepSz
  }
  
  # Start a half-step in from zlim[0]
  swatches <- zlim[1] + (stepSz/2)
  for(i in 2:length(palette)) {
    swatches[i] <- swatches[i-1] + stepSz
  }
  
  if(filename != "") {
    bitmap(filename, width=width, height=height, res=300)
    par(mar=c(2.5,1.2,0.5,1.2), mgp=c(1.5,0.5,0), lwd=0.5, cex=cex)
  }
  
  image(
    x=gridlines,
    z=as.matrix(swatches, ncol=1),
    zlim=zlim,
    col=palette,
    xlab=xlab,
    yaxt="n",
    xaxp=c(axis.limits, axis.intervals)
  )
  
  if(filename != "") {
    dev.off()
  }
}

cat("Drawing palette legend to:", my.args$PALFILE, "\n")
drawHeatmapScale(zlim=c(-1,1), palette=heatmap.pal, xlab="Correlation Coefficient", filename=my.args$PALFILE)

# Compute correlations within each cluster, and compile total list of within-cluster K-K, K-N, and N-N correlations
# Only do this if we have known and novel genes in this set
mmc.table <- mmc.table[row.names(crop.expr.cor),]
cat("Computing within-module correlation distributions.\n")
known.known.cor <- c()
novel.novel.cor <- c()
known.novel.cor <- c()
for(k in unique(mmc.table$Module)) {
  k.genes <- row.names(mmc.table)[mmc.table$Module==k]
  within.k.cor <- crop.expr.cor[k.genes,k.genes]
  # Separate known and novel gene IDs
  k.known <- grep("^FBgn",k.genes,value=T)
  k.novel <- grep("^XLOC",k.genes,value=T)
  if(length(k.known) > 1) {
    k.known.cor <- lapply(1:(length(k.known)-1), function(i){
      within.k.cor[k.known[i],k.known[(i+1):length(k.known)]]
    })
    known.known.cor <- c(known.known.cor, unlist(k.known.cor))
  }
  if(length(k.novel) > 1) {
    k.novel.cor <- lapply(1:(length(k.novel)-1), function(i){
      within.k.cor[k.novel[i],k.novel[(i+1):length(k.novel)]]
    })
    novel.novel.cor <- c(novel.novel.cor, unlist(k.novel.cor))
  }
  if((length(k.novel) > 0) & (length(k.known) > 0)) {
    known.novel.cor <- c(known.novel.cor, c(within.k.cor[k.known,k.novel], recursive=T))
  }
}
  
if(length(known.known.cor) > 0) {
  known.known.dens <- density(known.known.cor)
} else {
  known.known.dens <- list(x=c(), y=c())
}
if(length(novel.novel.cor) > 0) {
  novel.novel.dens <- density(novel.novel.cor)
} else {
  novel.novel.dens <- list(x=c(), y=c())
}
if(length(known.novel.cor) > 0) {
  known.novel.dens <- density(known.novel.cor)
} else {
  known.novel.dens <- list(x=c(), y=c())
}
y.max <- max(c(known.known.dens$y, novel.novel.dens$y, known.novel.dens$y))*1.05

cat("Drawing within cluster expression correlations to:", my.args$CORFILE, "\n")
pdf(my.args$CORFILE, width=5, height=3, pointsize=16)
par(mgp=c(1.7,0.6,0), mar=c(2.8,2.8,0.5,0.5), tcl=-0.4)
plot(x=0, y=0, type="n", xlim=c(-1,1), ylim=c(0,y.max), xlab="Within-Cluster Expression Correlation", ylab="Density")
lines(known.known.dens$x, known.known.dens$y, lwd=4, col="green3")
lines(novel.novel.dens$x, novel.novel.dens$y, lwd=4, col="cyan4")
lines(known.novel.dens$x, known.novel.dens$y, lwd=4, col="orange3")
dev.off()

