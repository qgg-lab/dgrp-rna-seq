#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 2/23/17
#
# expression_PCA.R
#
# Goal: Perform PCA on line means, output the PCs and what % of variance they explain for each gene
# 
# Usage:
# Rscript expression_PCA.R EXPR=subdir/data_line_means.txt [OUTPUT=subdir/output_prefix]
#  EXPR=      Path to file of line means
#  OUTPUT=    Path and prefix for all output files
#             Defaults to path of EXPR file (dropping genVar/, adding pca/) and taking the EXPR file name minus _line_means.txt
#  CORES=     Number of threads to use for parallelization
# 

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

options(stringsAsFactors=F)

# For parallelization of the linear regression steps
library(doMC)

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript expression_PCA.R EXPR=subdir/data_line_means.txt [OUTPUT=subdir/output_prefix]"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_F_line_means.txt")
# my.args <- c("EXPR=microbiome_species_filtered/genVar/combined_microbe_filtered_rpm_WolAdj_SexAvg_line_means.txt")

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

# Determine default OUTPUT value
default.output <- dirname(my.args$EXPR)
default.output <- sub("genVar$", "", default.output)
if(default.output == "") {
  default.output <- "."
}
if(!grepl("/$", default.output)) {
  default.output <- paste0(default.output, "/")
}
default.output <- paste0(default.output, "pca/", basename(my.args$EXPR))
default.output <- sub("[.]txt$", "", default.output)
default.output <- sub("[_]line[_]means$", "", default.output)

setDefault("OUTPUT", default.output)

# Get just the directory portion of OUTPUT, create if it doesn't exist
output.dir <- dirname(my.args$OUTPUT)
if(!file.exists(output.dir)) {
  dir.create(path=output.dir, recursive=T, showWarnings = F)
}

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


# --- LOAD INPUT FILES --- #

# Load the line means
cat("Loading expression values from:", my.args$EXPR, "\n")
expr.table <- read.table(my.args$EXPR, header=T, sep="\t", row.names=1, as.is=T)
colnames(expr.table) <- sub("^X","", colnames(expr.table))
cat("Loaded expression data for", nrow(expr.table), "features across", ncol(expr.table), "samples.\n\n")

# TO DO - Filter on FLAG columns here?

# Put flags in seperate data structure
if("FLAG" %in% colnames(expr.table)) {
  gene.flags <- expr.table$FLAG
  names(gene.flags) <- row.names(expr.table)
  expr.table <- expr.table[,setdiff(colnames(expr.table),"FLAG")]
} else {
  gene.flags <- rep("OK", times=nrow(expr.table))
  names(gene.flags) <- row.names(expr.table)
}

# Transpose the expression table and scale (NOTE: This automatically converts to matrix)
# NOTE: This means loading factors are meant to be applied to centered/scaled data, not the raw line means
expr.table <- t(expr.table)
expr.table <- scale(expr.table)
expr.center <- attr(expr.table, "scaled:center")
expr.sd <- attr(expr.table, "scaled:scale")


# --- Principle Component Analysis --- #

# Compute the variance
cat("Performing PCA analysis...\n")
expr.prcomp <- prcomp(expr.table, retx=T)
pc.var <- expr.prcomp$sdev ** 2
pc.pvar <- pc.var / sum(pc.var)
pc.cvar <- unlist(lapply(1:length(pc.pvar), function(i){sum(pc.pvar[1:i])}))
pc.table <- expr.prcomp$x

pc.summary <- cbind(data.frame(PC=colnames(pc.table), Var=pc.var, PropVar=pc.pvar, CumPropVar=pc.cvar))
cat("Top PCs:\n")
print(head(pc.summary))

# Limit to PCs that explain at least 0.5% of variance for now
pc.summary <- pc.summary[pc.summary$PropVar >= 0.005,]
cat("Keeping", nrow(pc.summary), "PCs individually explaining at least 0.5% of variance,\n")
cat("Cumulative % variance explained =", round(sum(pc.summary$PropVar)*100), "%\n")
pc.table <- pc.table[,pc.summary$PC]

# Build the general LM formula
# (This is necessary to get ANOVA to give individual P-value for each PC)
expr.pc.formula <- as.formula(paste("EXPR ~", paste(pc.summary$PC, collapse=" + ")))


# Loop over all genes, fit lm, do anova, compute % var explained
cat("Testing % Var Explained and ANOVA P-values for",length(pc.summary$PC),"PCs against", nrow(expr.table), "features.\n")
expr.pc.results <- foreach(gene=colnames(expr.table), .combine=rbind) %dopar% {
  # Combine this gene expression profile with these PCs
  gene.lm.input <- as.data.frame(cbind(EXPR=expr.table[,gene], pc.table))
  # Fit Linear Model over just the top PCs
  gene.pc.lm <- lm(expr.pc.formula, data=gene.lm.input)
  gene.pc.coeff <- gene.pc.lm$coefficients[pc.summary$PC]
  
  # Perform ANOVA to get P-value for each PC against this gene expression profile
  gene.pc.anova <- as.data.frame(anova(gene.pc.lm))
  
  # Compute % variance explained by each PC
  gene.var <- var(expr.table[,gene])
  gene.pc.propvar <- unlist(lapply(pc.summary$PC, function(pc){
    var.wout.pc <- var(expr.table[,gene] - (gene.pc.coeff[pc] * pc.table[,pc]))
    (gene.var - var.wout.pc)/gene.var
  }))
  names(gene.pc.propvar) <- pc.summary$PC
  gene.pc.propvar[gene.pc.propvar < 0] <- 0
  
  # Interpolate PropVar and PVal for each PC
  unlist(lapply(pc.summary$PC, function(pc){
    pc.values <- c(gene.pc.propvar[pc], gene.pc.anova[pc,"Pr(>F)"])
    names(pc.values) <- paste(pc, c("Perc.Var","PVal"), sep=".")
    return(pc.values)
  }))
}
row.names(expr.pc.results) <- colnames(expr.table)
expr.pc.results <- as.data.frame(expr.pc.results)

# Compute FDR adjustment for each PC PVal column
for(pval in grep("[.]PVal$", colnames(expr.pc.results), value=T)) {
  expr.pc.results[,paste0(pval,".FDR")] <- p.adjust(expr.pc.results[,pval], method = "BH")
}

# Insert the loading variable columns into this table as well
PC.loading <- expr.prcomp$rotation[,pc.summary$PC]
stopifnot(nrow(expr.pc.results)==nrow(PC.loading))
stopifnot(all(row.names(expr.pc.results)==row.names(PC.loading)))
colnames(PC.loading) <- paste0(colnames(PC.loading), ".Loading")
expr.pc.results <- cbind(expr.pc.results, PC.loading)

# Reorder columns
col.order <- paste(rep(pc.summary$PC, each=4), rep(c("Loading","Perc.Var","PVal","PVal.FDR"), times=length(pc.summary$PC)), sep=".")
stopifnot(all(col.order %in% colnames(expr.pc.results)))
expr.pc.results <- expr.pc.results[,col.order]

# For each PC, compute the number of genes with significant ANOVA P-value (5% FDR)
# And the mean % Var Explained
pc.summary[,"SigGenes"] <- unlist(lapply(pc.summary$PC, function(pc){
  sum(expr.pc.results[,paste0(pc,".PVal.FDR")] <= 0.05)
}))
pc.summary[,"Avg.Gene.PropVar"] <- unlist(lapply(pc.summary$PC, function(pc){
  mean(expr.pc.results[,paste0(pc,".Perc.Var")])
}))
# Avg.Gene.PropVar should always be very close to original PropVar reported by PCA
if(any(abs(pc.summary$PropVar - pc.summary$Avg.Gene.PropVar) > 0.001)) {
  cat("WARNING: Not all Avg.Gene.PropVar values are = original PCA PropVar.\n")
} else {
  # Drop the column if it's redundant
  pc.summary <- pc.summary[,setdiff(colnames(pc.summary),"Avg.Gene.PropVar")]
}

# --- Draw Summary Plots --- #
# Plot PropVar, CumPropVar, and SigGenes for these

propvar.plot <- paste0(my.args$OUTPUT, "_PropVar_barplot.pdf")
cat("Drawing % Variance Explained plot to:", propvar.plot, "\n")
pdf(propvar.plot)
barplot(pc.summary$PropVar*100, col="lightblue", names=1:nrow(pc.summary), xlab="Principal Components", ylab="% of Variance Explained")
dev.off()

cumvar.plot <- paste0(my.args$OUTPUT, "_CumPropVar_barplot.pdf")
cat("Drawing Cumulative % Variance Explained plot to:", cumvar.plot, "\n")
pdf(cumvar.plot)
barplot(pc.summary$CumPropVar*100, col="indianred2", names=1:nrow(pc.summary), xlab="Principal Components", ylab="Cumulative % Variance Explained")
dev.off()

siggenes.plot <- paste0(my.args$OUTPUT, "_SigGenes_barplot.pdf")
cat("Drawing Significant Gene plot to:", siggenes.plot, "\n")
pdf(siggenes.plot)
barplot(pc.summary$SigGenes, col="lightgreen", names=1:nrow(pc.summary), xlab="Principal Components", ylab="Number of Significant Genes")
dev.off()

combined.plot <- paste0(my.args$OUTPUT, "_PropVar_SigGenes_plot.pdf")
cat("Drawing combined plot to:", combined.plot, "\n")
pdf(combined.plot)
propvar.range <- c(0,ceiling(max(pc.summary$PropVar*100)/5)*5)
propvar.ticks <- ((propvar.range[1]/5):(propvar.range[2]/5))*5
propvar.col <- "blue3"
par(mar=c(5,5,2,5))
plot(x=1:nrow(pc.summary), y=pc.summary$PropVar*100, type="o", pch=19, col=propvar.col, bty="n", xlab="Principal Components", xaxp=c(1, nrow(pc.summary), nrow(pc.summary)-1), ylim=propvar.range, yaxt="n", ylab="")
axis(side=2, col=propvar.col, col.axis=propvar.col, at=propvar.ticks)
mtext("% Variance Explained", side=2, line=3, col=propvar.col)
par(new=T)
siggenes.range <- c(0,ceiling(max(pc.summary$SigGenes)/2000)*2000)
siggenes.ticks <- ((siggenes.range[1]/2000):(siggenes.range[2]/2000))*2000
siggenes.col <- "forestgreen"
plot(x=1:nrow(pc.summary), y=pc.summary$SigGenes, type="o", pch=19, col=siggenes.col, bty="n", xlab="", xaxt="n", ylim=siggenes.range, yaxt="n", ylab="")
axis(side=4, col=siggenes.col, col.axis=siggenes.col, at=siggenes.ticks)
mtext("# of Significant Genes", side=4, line=3, col=siggenes.col)
dev.off()


# --- Write Tables --- #

# Output the table of PC values over each line (transpose back to match line means table format)
pc.table <- as.data.frame(t(pc.table))
pc.table <- cbind(GENE=row.names(pc.table), pc.table)
pc.table.file <- paste0(my.args$OUTPUT, "_PCs.txt")
cat("Writing PC vectors to:", pc.table.file, "\n")
write.table(pc.table, pc.table.file, row.names=F, sep="\t", quote=F)

# Output the table of PC summary stats
pc.summary.file <- paste0(my.args$OUTPUT, "_PC_stats.txt")
cat("Writing PC summary statistics to:", pc.summary.file, "\n")
write.table(pc.summary, pc.summary.file, row.names=F, sep="\t", quote=F)

# Output the tabel of PC results for each gene (make sure row names = genes)
expr.pc.results <- cbind(GENE=row.names(expr.pc.results), expr.pc.results)
expr.pc.file <- paste0(my.args$OUTPUT, "_PC_gene_results.txt")
cat("Writing Gene x PC model results to:", expr.pc.file, "\n")
write.table(expr.pc.results, expr.pc.file, row.names=F, sep="\t", quote=F)
