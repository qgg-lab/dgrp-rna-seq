#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 2/7/17
#
# This script takes filtered line means and runs WGCNA on all genes to get modules and network structure
#
# Usage:
# Rscript run_wgcna.R [OPTIONS] EXPR=line_means_wgcna.txt
#  EXPR=  The file containing line mean expression values to filter and prep for downstream analysis
#  CORES= The number of CPUs to use for parallelization (defaults to 1)
#   NOTE: If submitting your job through SLURM, make sure to set -c option as well
#
# Output will be in same directory as file specified by EXPR, with similar names (everything but the _line_means_wgcna.txt part)
#

# TO DO: Add command line params to adjust the soft-power, and the correlation method
# (There should be a way to tag the output files as separate from defaults as well)
#

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

# TO DO: PUT URL TO INSTALL INSTRUCTIONS HERE
# Library is already installed on Hyperion under ~ljeveret/Tools/R-3.1.1/bin/R
library(WGCNA)

options(stringsAsFactors=F)

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript run_wgcna.R [OPTIONS] EXPR=line_means_wgcna.txt"
my.args <- commandArgs(trailingOnly=T)
# TEMP TESTING:
# my.args <- c("EXPR=known_all_novel_genes/wgcna/combined_samples_known_novel_fpkm_VR_WolAdj_M_line_means_wgcna.txt", "CORES=4")
# my.args <- c("EXPR=microbiome_species_filtered/wgcna/combined_microbe_filtered_rpm_WolAdj_SexAvg_line_means_wgcna.txt")

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

# EXPR=  The file containing line mean expression values to filter and prep for downstream analysis
# Make sure EXPR param exists
if(!("EXPR" %in% names(my.args))) {
  stop("Missing EXPR parameter. ", usageStr)
}

# Determine output directory and stub
my.args$OUTDIR <- paste0(dirname(my.args$EXPR), "/")
my.args$OUTSTUB <- basename(my.args$EXPR)
my.args$OUTSTUB <- sub("[.]txt$", "", my.args$OUTSTUB)
my.args$OUTSTUB <- sub("[_]wgcna$", "", my.args$OUTSTUB)
my.args$OUTSTUB <- sub("[_]line[_]means$", "", my.args$OUTSTUB)

setDefault("CORES",as.numeric(NA))
if(is.na(my.args$CORES)) {
  # If NA, see if environment variable SLURM_JOB_CPUS_PER_NODE exists
  slurm.cpus <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE")
  if(slurm.cpus != "") {
    slurm.cpus <- as.numeric(slurm.cpus)
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
if(is.character(my.args$CORES) & (my.args$CORES == "ALL")) {
  my.args$CORES <- as.numeric(detectCores()-1)
  cat("Setting CORES =", my.args$CORES, "\n")
}
# Otherwise, attempt to convert to integer
if(!is.numeric(my.args$CORES)) {
  my.args$CORES <- as.numeric(my.args$CORES)
  cat("Setting CORES =", my.args$CORES, "\n")
}

# Enable multi-threading (TO DO: make this adjustable)
if(my.args$CORES > 1) {
  enableWGCNAThreads(nThreads=my.args$CORES)
}

cat("\n")


# Load the line means (should already be in sample by gene format, tab-delimited)
cat("Loading line means from:", my.args$EXPR, "\n")
lineMeans <- read.table(my.args$EXPR, header=T, sep="\t", row.names=1)
cat("Loaded line means for", ncol(lineMeans), "features across", nrow(lineMeans), "lines.\n\n")


# TO DO: Can generate my own similarity matrix (just needs to have values in range [0,1])
#        and use pickSoftThreshold.fromSimilarity and adjacency.fromSimilarity functions 
#        to determine optimal power and do the transformation
#        
# One idea I had was to rank-transform the correlations to normalize connectivity across the network
# e.g., similarity(i,j) = (rank(cor(i,j),cor(i,*))+rank(cor(i,j),cor(j,*)))/2
#   Where rank(x,X) is the normalized rank of x in vector X (returns 0 if x = min(X), returns 1 if x = max(X))
#   So this would pare down some of the connectivity in very densely correlated parts of the network and boost connectivity in less correlated parts
#   This might actually work quite well with TOM - it says that in heavily connected parts of the network, 
#   two nodes need to have their strongest connections to the SAME nodes in order to end up in the same module
#
# Another idea discussed in Zhang and Horvath 2005 is to use the significance of the correlation
#  which can be obtained using the Fisher transformation (Davidson, et al. 2001)
#  or using permutation tests (Butte and Kohane, 2000 and Carter et al. 2004)
# However, they only present this idea for getting hard-threshold similarity matrix, e.g. picking a significance cut-off
# So would need some way to convert this to a soft-threshold similarity matrix (1-pval?)
# OR could do (1-pval)*cor to penalize noisier correlation coefficents?
#


# --- Determine soft-threshold --- #

# This is the range suggested in the first tutorial, may need to adjust this?
powers <- c(c(1:10), seq(from=12, to=20, by=2))
cat("Picking soft-threshold power, testing selected powers from", min(powers), "to", max(powers), "\n")

# TESTING DIFFERENT CORRELATION METHODS HERE
# USES PEARSON BY DEFAULT, BUT SPEARMAN MIGHT BE BETTER (BUT IS MUCH SLOWER)
# Note that their built-in optimizations are only available for Pearson Correlation, so Spearman is much slower
# ALSO: Spearman might work better if all gene expression profiles are centered and maybe scaled first?
# Could try converting expr matrix to rank matrix first, and then use Pearson Correlation here?
# The 2005 paper (Zhang and Horvath) also suggests using jackknifing to make the correlation robust to outliers but there is no implementation of this in their functions
# However, could do that externally and then use the pickSoftThreshold.fromSimilarity function to provide directly the matrix of jackknifed correlations
# Jackknife is similar to bootstrap except that instead of sampling with replacement, you randomly drop 1 or more samples on each subsampling
# More info here: https://en.wikipedia.org/wiki/Resampling_(statistics)#Jackknife
# TO DO: Should also setting moreNetworkConcepts=T to get more network params computed by this function
# TO DO: Look at sft$powerEstimate for the automatically selected value of optimal power here!
# TO DO: Also plot the truncated version of scale-free fit?
# Some notes based on reading Zhang and Horvath 2005:
#  The scale-free topology metric used here is to first plot log10(p(k)) vs log10(k)
#   To do this, there is some binning of k values (since they're not actually integers using soft thresholding)
#  Then then try to fit a straight line to these values, and the R^2 of that line is the metric (closer to 1 = better fit)
#  They also fit a more complicated model: p(k) ~ k^(-y)*exp(-ak) that has two fit params, y and a, which produces the truncated R^2
#  Also, the slope of the best-fit line is either negative (expected) or positive (weird), 
#   which is why they multiply by the opposite sign of the slope.
#   In effect, this invalidates good R^2 values when the slope goes the wrong way.
#  They generally think that an R^2 >= 0.85 (or 0.8?) is sufficiently scale-free, so long as the slope is negative.
sft <- pickSoftThreshold(lineMeans, powerVector=powers, verbose=5, corOptions=list(use='p', method='pearson'))

# The scale independence looks pretty wonky for Female line means
# Will try to use mean connectivity instead

# DEPRECATING THIS PART...
# Try fitting an exponential regression on mean.k.
# mean.k.exp.reg <- lm(log(mean.k.) ~ Power, data=sft$fitIndices)
#
# Make a prediction function
# mean.k.predict <- function(Power) {
#  exp((Power * mean.k.exp.reg$coefficients["Power"]) + mean.k.exp.reg$coefficients["(Intercept)"])
# }
# Make a first derivative (slope) function
# mean.k.slope <- function(Power) {
#  mean.k.exp.reg$coefficients["Power"] * exp((Power * mean.k.exp.reg$coefficients["Power"]) + mean.k.exp.reg$coefficients["(Intercept)"])
# }
# Make a second derivative (inflection) function
# mean.k.inflect <- function(Power) {
#  (mean.k.exp.reg$coefficients["Power"]^2) * exp((Power * mean.k.exp.reg$coefficients["Power"]) + mean.k.exp.reg$coefficients["(Intercept)"])
# }
# mean.k.predict.x <- (1:200)/10
# mean.k.predict.y <- mean.k.predict(mean.k.predict.x)
# mean.k.slope.y <- mean.k.slope(mean.k.predict.x)
# mean.k.inflect.y <- mean.k.inflect(mean.k.predict.x)
#
# For each mean-connectivity point after power=1, compute the slope from previous power
# sft$fitIndices[1,"mean.k.slope"] <- as.numeric(NA)
# for(i in 2:nrow(sft$fitIndices)) {
#  sft$fitIndices[i,"mean.k.slope"] <- (sft$fitIndices[i,"mean.k."] - sft$fitIndices[i-1,"mean.k."]) / (sft$fitIndices[i,"Power"] - sft$fitIndices[i-1,"Power"])
# }

# TO DO: When plotting mean.k. it would make sense to scale this to the number of nodes in the network?
#       i.e. range would be 0% (no connectivity) to 100% (full connectivity)
# TO DO: Other criteria for selecting power (suggested in Zhang and Horvath 2005) are:
#       Try to maximize connectivity so that strong modules can be found (within the bounds of other criteria)
#       The slope of the regression line for the scale-free model should be closer to -1 (this is 3rd column in fitIndices)
#       

# Use the WGCNA automatic power threshold if available
if(!is.na(sft$powerEstimate)) {
  use.power <- sft$powerEstimate
  cat("Using soft-power =", use.power, "determined by WGCNA\n")
} else {
  # This was my original approach
  # Compute the range of mean-connectivity values, pick a threshold that is within 5% (of the range) of the lower value
  # mean.k.range <- diff(range(sft$fitIndices[,"mean.k."]))
  # mean.k.threshold <- min(sft$fitIndices[,"mean.k."]) + (0.05*mean.k.range)
  # Soft-power threshold is the lowest power that gets us below mean.k.threshold
  # use.power <- min(sft$fitIndices[sft$fitIndices[,"mean.k."] <= mean.k.threshold,"Power"])
  # cat("Using soft-power threshold =", use.power, "based on mean connectivity threshold <=", mean.k.threshold, "\n")
  
  # Pick lowest power to get within 10% of max signed SFT R2
  sft.cutoff <- max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]) * 0.9
  use.power <- min(sft$fitIndices[which((-sign(sft$fitIndices[,3])*sft$fitIndices[,2]) >= sft.cutoff),1])
  cat("Using soft-power =", use.power, "based on Scale-Free Topology R^2 curve.\n")
}

# Draw the PDF graph
pdf.file <- paste0(my.args$OUTDIR, my.args$OUTSTUB, "_wgcna_powers.pdf")
cat("Plotting Scale-Free Topology statistics to:", pdf.file, "\n")
pdf(pdf.file, width=9, height=5)
par(mfrow=c(1,2))
cex1=0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
points(sft$fitIndices[use.power,1], -sign(sft$fitIndices[use.power,3])*sft$fitIndices[use.power,2], pch=1, cex=2, col="purple")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# abline(h=mean.k.threshold, col="blue")
points(sft$fitIndices[use.power,1], sft$fitIndices[use.power,5], pch=1, cex=2, col="purple")
# DEPRECATED: Was looking at exponential regression and derivatives of best-fit curve but that turned out not to be so useful
# lines(x=mean.k.predict.x, y=mean.k.predict.y, col="blue")
# lines(x=mean.k.predict.x, y=mean.k.slope.y, col="green")
#lines(x=mean.k.predict.x, y=mean.k.inflect.y, col="purple")
dev.off()
cat("\n")


# --- Run the network construction --- #

# Manual step-by-step approach gives more parameter flexibility
# It looks like there is no blocksize parameter in any of these functions, so that is only an issue with the automatic construction function?
# The original 2005 Zhang and Horvath paper also suggests using the sigmoid function for converting similarity matrix to adjacency matrix
#   But there is no implementation of that here, and paper suggests there is no advantage compared to power function
# TO DO: Should just try using power =1 here, i.e. adjacency = similarity
cat("Computing adjacency matrix for all genes based on line means.\n")
adjacency = adjacency(lineMeans, power = use.power, corOptions="use = 'p', method = 'pearson'")

# Turn adjacency into topological overlap
cat("Computing Topological Overlap Map.\n")
TOM = TOMsimilarity(adjacency, verbose=3)
dissTOM = 1-TOM

# Call the hierarchical clustering function
cat("Computing hclust tree from TOM distance matrix.\n")
geneTree = hclust(as.dist(dissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
# Default values were minModuleSize=30, deepSplit=2
# Increasing deepSplit to 4 (max?) and decreasing minModuleSize to 20 puts a few more unassigned genes in modules
# and splits up some of the modules, but doesn't really do anything about the primary large module of ~5k genes
# Worth taking a close look at this function documentation, there is a lot of fine control params available

# Appropriate minModuleSize depends on how many things are being clustered
# TO DO: Make this user adjustable?
if(ncol(lineMeans) >= 200) {
  minModuleSize <- 20
} else {
  minModuleSize <- 4
}

# Module identification using dynamic tree cut:
cat("Searching for modules of >=", minModuleSize, "features.\n")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize, verbose = 3);

cat("Found", length(setdiff(unique(dynamicMods), "0")), "clusters with sizes (0=Unclustered):\n")
table(dynamicMods)
cat("\n")

# Plot the dendrogram and module colors
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

# Plot the dendrogram and colors underneath
# SKIP - GETS DONE AGAIN AFTER MERGING STEP
# pdf.file <- paste0(my.args$OUTDIR, my.args$OUTSTUB, "_wgcna_dendro_modules.pdf")
# cat("Drawing dendrogram with module color map to:",pdf.file,"\n")
# pdf(pdf.file, width=12, height=9)
# plotDendroAndColors(geneTree, dynamicColors, "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Gene dendrogram and module colors")
# dev.off()
# cat("\n")

# Calculate eigengenes
cat("Calculating and drawing eigengenes.\n")
MEList = moduleEigengenes(lineMeans, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Make this adjustable or skip the merge step?
MEDissThres = 0.25

pdf.file <- paste0(my.args$OUTDIR, my.args$OUTSTUB, "_wgcna_eigen_dendro.pdf")
cat("Drawing clustering of module eigengenes to:",pdf.file,"\n")
pdf(pdf.file, width=7, height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()
cat("\n")

# Call an automatic merging function
# TO DO: Should set to use Spearman correlation here too?  It's not clear to me if this correlation is between eigengenes or between individual genes?
merge = mergeCloseModules(lineMeans, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Re-order the cluster numbers starting from 1, with biggest first (similar to what cutreeDynamic does)
colorSz <- table(mergedColors)
colorSz <- colorSz[setdiff(names(colorSz),"grey")]
colorSz <- sort(colorSz, decreasing=T)
colorID <- 0:(length(colorSz))
names(colorID) <- c("grey",names(colorSz))
mergedMods <- colorID[mergedColors]
stopifnot(sum(is.na(mergedMods))==0)
cat("Merged into", length(setdiff(unique(mergedMods), "0")), "clusters with sizes (0=Unclustered):\n")
table(mergedMods)
cat("\n")
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# Draw dendrogram again with both dynamic and merged cluster colors
pdf.file <- paste0(my.args$OUTDIR, my.args$OUTSTUB, "_wgcna_dendro_modules.pdf")
cat("Drawing dendrogram with merged module color map to:",pdf.file,"\n")
pdf(pdf.file, width=12, height=9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
cat("\n")

# Draw a PNG file using random subset of genes
# (Trying N=1000)
# nSelect <- 1000
# For reproducibility, set random seed
# set.seed(8210678)
# select <- sample(ncol(lineMeans), size=nSelect)
# selectTOM <- dissTOM[select,select]
# selectTree <- hclust(as.dist(selectTOM), method="average")
# selectColors <- dynamicColors[select]
# Open a graphical window
# DEPRECATED - NO BITMAP DRAWING ON CLUSTER, SUPER ANNOYING!!!!!!
# TO DO: PUT THIS IN A SEPARATE PLOTTING SCRIPT
# png.file <- paste0(my.args$OUTDIR, my.args$OUTSTUB, "_wgcna_tom_heatmap_selected.png")
# cat("Drawing TOM heatmap for selected subset of genes to:", png.file, "\n")
# bitmap(png.file, width=9, height=9)
# # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# # the color palette; setting the diagonal to NA also improves the clarity of the plot
# plotDiss = selectTOM^7;
# diag(plotDiss) = NA;
# TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
# dev.off()
# cat("\n")

# Try drawing map of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(lineMeans, dynamicColors)$eigengenes
MET = orderMEs(MEs)
# Plot the relationships among the eigengenes and the trait
pdf.file <- paste0(my.args$OUTDIR, my.args$OUTSTUB, "_wgcna_eigengenes.pdf")
cat("Drawing Eigengene heatmap to:", pdf.file, "\n")
pdf(pdf.file, width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()

# Create a table with genes and module assignments
module.table <- data.frame(Gene=colnames(lineMeans), Module=as.integer(mergedMods))
module.table[,"Entry Index"] <- 1:ncol(lineMeans)
# For each numbered module (not 0), compute the Degree and Average Degree (use the TOM adjacency matrix)
module.table[,"Average Degree"] <- as.numeric(NA)
module.table[,"Degree"] <- unlist(lapply(1:nrow(module.table), function(i){
  i.mod <- module.table[i,"Module"]
  i.mod.other <- setdiff(module.table[module.table$Module==i.mod,"Entry Index"], i)
  mean(TOM[i,i.mod.other])
}))
# Set Degree = 0 for unclustered genes
module.table[module.table$Module==0,"Degree"] <- 0
# Sort by cluster, degree
module.table <- module.table[order(module.table$Degree, decreasing=T),]
module.table <- module.table[order(module.table$Module),]
# Compute Avg Degree for each cluster, sort within each cluster by degree
for(k in unique(module.table$Module)) {
  module.table[module.table$Module==k,"Average Degree"] <- mean(module.table[module.table$Module==k,"Degree"])
}
# Re-sort by Average Degree (this will put cluster 0 at the end)
module.table <- module.table[order(module.table[,"Average Degree"], decreasing=T),]
# Re-number clusters based on Avg Degree order
new.cluster.order <- setdiff(unique(module.table$Module), 0)
new.cluster.map <- 1:length(new.cluster.order)
names(new.cluster.map) <- as.character(new.cluster.order)
new.cluster.map["0"] <- 0
module.table$Module <- new.cluster.map[as.character(module.table$Module)]

# Output this table in CSV format so that it can be handled the same way as MMC output
output.file <- paste0(my.args$OUTDIR, my.args$OUTSTUB, "_wgcna.csv")
cat("Writing module assignments to:", output.file, "\n")
write.table(module.table, output.file, sep=",", quote=F, row.names=F)

# TO DO: Also check out the propVarExplained function to see what proportion of each gene's expression variance is explained by corresponding eigengene

# TO DO: May also be worth looking at cluster coefficients here to look at overall modularity
# (see discussion at end of Zhang and Horvath 2005)

# TO DO: Could load phenotype data here and look for eigengenes that are predictive
# Rather than using simple correlation, Zhang and Horvath (2005) use univariate regression models 
# (Cox only applies to survival data, so not approporiate here, but simple lm/anova should work?)
# Worth considering this as replacement for the simple correlation analysis I'm currently using
# NOTE: They use -log10(p) as final metric, where p denotes the univariate regression p-value


cat("\n\nScript completed successfully!\n")
print(proc.time())
cat("\n")
quit("no")


# --- OLD CODE --- #

# THIS WAS THE AUTOMATED VERSION, USING STEPWISE (ABOVE) INSTEAD TO HAVE MORE CONTROL...
# TO DO: May need to auto-adjust maxBlockSize for more genes, and/or make this user configurable
net = blockwiseModules(lineMeans, power=use.power, maxBlockSize=10000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "combined_samples_known_novel_fpkm_VR_WolAdj_F_TOM", 
                       verbose = 3)

# How many modules, and what are their sizes?
table(net$colors)
# Ick, this is still terrible!

# TO DO: Extract the first PC from lineMeans and treat it as a trait, see how it correlates with overall module structure here!

# Try drawing the dendrogram
pdf("combined_samples_known_novel_fpkm_VR_WolAdj_F_wgcna_dendro_modules.pdf", width=12, height=9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Try drawing the topological overlap heatmap
# This is more akin to the MMC correlation heatmap

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(lineMeans, power = use.power);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

# Draw a PNG file - this FAILS due to memory limits (may work on exclusive node if R run w/ params to increase mem limit?)
# bitmap("combined_samples_known_novel_fpkm_VR_WolAdj_F_wgcna_tom_heatmap_all.png", width=9, height=9)
# TOMplot(plotTOM, net$dendrograms[[1]], labels2colors(net$colors), main = "Network heatmap plot, all genes")
# dev.off()

# THIS CAN'T RUN ON CLUSTER - NEED TO PUT IN A SEPERATE SCRIPT OR MAKE IT PART OF THE DOWNSTREAM PLOTTING SCRIPT
# # Draw a PNG file using random subset of genes
# # (Trying N=1000)
# nSelect <- 1000
# # For reproducibility, set random seed
# set.seed(8210678)
# select <- sample(ncol(lineMeans), size=nSelect)
# selectTOM <- dissTOM[select,select]
# selectTree <- hclust(as.dist(selectTOM), method="average")
# selectColors <- labels2colors(net$colors)[select]
# # Open a graphical window
# bitmap("combined_samples_known_novel_fpkm_VR_WolAdj_F_wgcna_tom_heatmap_selected.png", width=9, height=9)
# # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# # the color palette; setting the diagonal to NA also improves the clarity of the plot
# plotDiss = selectTOM^7;
# diag(plotDiss) = NA;
# TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
# dev.off()

# Try drawing map of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(lineMeans, labels2colors(net$colors))$eigengenes
# See FemaleLiver-05-Visualization tutorial for adding in traits of interest here
MET = orderMEs(MEs)
# Plot the relationships among the eigengenes and the trait
pdf("combined_samples_known_novel_fpkm_VR_WolAdj_F_wgcna_eigengenes.pdf", width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()


