#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 9/1/16
#
# Script to take library count data and basic info and check for any samples
# with incorrect sex label or possible opposite sex contamination
#
# USAGE: validate_sex.R [OPTIONS] BATCH1 [BATCH2 ...]
#  --samples=FILE specifies the sample info file (Default: sample_info.txt)
#  --path=PATH specifies the top-level directory containing HTSeq-count results (Default: htseq/)
#  --suffix= specifies the part of count file names that comes after the library ID (Default: _STAR_counts.txt)
#  --qcpath=PATH specifies the directory to output QC plots (Default: sample_QC_figures/)
#
# The script loads the gene counts per library, 
# does a simple RPM normalization using total counts, 
# followed by log2 conversion, then 2 primary analyses:
# 1) PCA analysis assuming that PC1 will be sex
# 2) libraries are then split by sex and expr scaled w/in sex
#     The script then plots the z-score expr of different quantiles
#     Problematic samples and batches should stand out on this plot
#
# QC plots are also output to [QCPATH]/libsex_[STUB]_....pdf
#
# TO DO: Could try running on only the top 10k genes by mean expression to see if results are a bit more robust?
#

# Needed for fast stitching of HTSeq-count files into one big table:
library(foreach)
# Neded for linear discriminant analysis
library(MASS)
# EdgeR is not currently used, but if running into problems where PC1 does not capture main sex difference
# Then perhaps it would best to bring back the TMM normalization?
# library(edgeR)
options(stringsAsFactors=F)

# -- Fixed Parameters -- #
# These could be made adjustable in the command-line args
#
# Change pseudo-count here, this sets the minimum RPM values when converting to RLE
# TO DO: Make this adjustable or auto-determine an appropriate value for this!
pseudo.count <- 0.001


# --- PARSE COMMAND LINE PARAMS --- #

# Take the list of batches to analyze from command line params
usageStr <- "USAGE:\nRscript genotype_samples.R [OPTIONS] BATCH1 [BATCH2 ...]\n  --samples=FILE specifies the sample master table file (Default: sample_info.txt)\n  --path=PATH specifies the top-level directory containing allele counts (Default: allele_counts/)"
myArgs <- commandArgs(trailingOnly=T)
# Interactive code testing: 
# myArgs <- c("batch_160311A","batch_160311B")
# Separate out option args (begin with --) from standard args
myParams <- myArgs[grepl("^--",myArgs)]
batches <- myArgs[!grepl("^--",myArgs)]

if(length(batches) == 0) {
  # Can't run without at least one batch!
  stop("Must specify one or more batches to verify genotypes. ", usageStr)
}

# Remove trailing .txt from these files (just in case)
batches <- sub("[.]txt$", "", batches)

cat("Analyzing batches:", batches, "\n\n")

# Parse myParams
myParams <- strsplit(myParams, split="=", fixed=T)
names(myParams) <- unlist(lapply(myParams, function(x){sub("^--","",x[1])}))
myParams <- lapply(myParams, "[", 2)

# Fill in defaults
if(is.null(myParams$samples)) {
  cat("Setting --samples=sample_info.txt by default.\n")
  myParams$samples <- "sample_info.txt"
}

if(is.null(myParams$path)) {
  cat("Setting --path=htseq/ by default.\n")
  myParams$path <- "htseq/"
}
if(!grepl("/$",myParams$path)) {
  myParams$path <- paste0(myParams$path, "/")
}

if(is.null(myParams$suffix)) {
  cat("Setting --suffix=_STAR_counts.txt by default.\n")
  myParams$suffix <- "_STAR_counts.txt"
}

if(is.null(myParams$qcpath)) {
  cat("Setting --qcpath=sample_QC_figures/ by default.\n")
  myParams$qcpath <- "sample_QC_figures/"
}
# If QCPATH is missing trailing "/" add it on
if((!grepl("/$", myParams$qcpath)) & (myParams$qcpath != "")) {
  myParams$qcpath <- paste0(myParams$qcpath, "/")
}

cat("\n")




# --- Load Input Files --- #

# Load all batch files and append into a single table
library.table <- NULL
for(batch in batches) {
  batch.file <- paste0(batch, ".txt")
  if(!file.exists(batch.file)) {
    stop("Missing batch file: ", batch.file)
  }
  batch.table <- read.table(batch.file, sep="\t", row.names=1)
  colnames(batch.table) <- c("SAMPLE", "FLOWCELL", "LANE", "BARCODE", "DATE")
  batch.table[,"BATCH"] <- batch
  library.table <- rbind(library.table, batch.table)
}

# Load sample_info table
if(!file.exists(myParams$samples)) {
  stop("Missing sample info file: ", myParams$samples)
}
sample.info <- read.table(myParams$samples, header=T, sep="\t", row.names=1)

# Make sure all sample names in batch files are also in sample.info
if(!all(library.table$SAMPLE %in% row.names(sample.info))) {
  missing.samples <- setdiff(unique(library.table$SAMPLE), row.names(sample.info))
  stop(myParams$samples, " is missing the following sample IDs: ", paste(missing.samples, sep=" "))
}

# Make sure sample info has SEX column
if(!("SEX" %in% colnames(sample.info))) {
  stop(myParams$samples, " does not contain SEX column to validate.")  
}

# Map SEX column onto library.table
library.table[,"SEX"] <- sample.info[library.table$SAMPLE,"SEX"]

# Create list of count files for all libraries
library.files <- paste0(myParams$path, library.table$BATCH, "/", row.names(library.table), myParams$suffix)
names(library.files) <- row.names(library.table)

# Make sure all library files exist
libfile.exists <- unlist(lapply(library.files, file.exists))
if(!all(libfile.exists)) {
  stop("Missing library count files: ", paste(library.files[!libfile.exists], collapse=" "))
}

# Load the libraries
lib.counts <- foreach(lib=library.files, .combine=cbind) %do% read.table(lib, sep="\t", row.names=1)
colnames(lib.counts) <- names(library.files)

# Drop the "__*" rows from HTseq
lib.counts <- lib.counts[!grepl("^__", row.names(lib.counts)),]
cat("Loaded raw read counts for", nrow(lib.counts), "gene features across", ncol(lib.counts), "libraries.\n\n")

# --- Normalization and Log2 Transform --- #

cat("Performing RPM normalization and conversion to Log2 scale.\n")

# Normalize to RPM based on library gene feature count totals
lib.count.sums <- apply(lib.counts, 2, sum)
rpm.table <- lib.counts
for(j in 1:ncol(lib.counts)) {
  rpm.table[,j] <- lib.counts[,j] * (10**6) / lib.count.sums[j]
}

# Compute average RPM for each gene across libraries
gene.avg.rpm <- apply(rpm.table, 1, mean)
rpm.table <- rpm.table[gene.avg.rpm >= 1,]
cat("Using", nrow(rpm.table), "genes w/ avg RPM >= 1 for further validation.\n\n")

# Make sure libraries are in the same order in both tables
stopifnot(sum(is.na(library.table))==0)
stopifnot(sum(is.na(rpm.table))==0)
stopifnot(all(row.names(library.table)==colnames(rpm.table)))

female.libs <- row.names(library.table)[library.table$SEX=="F"]
male.libs <- row.names(library.table)[library.table$SEX=="M"]

cat(length(female.libs), "Libraries are labeled as Female and", length(male.libs), "Libraries are labeled as Male.\n")

cat("\n")

# DEPRECATED - TMM NORMALIZATION
# # Normalize to RPM based on lib count sums, but also use EdgeR TMM lib size normalization
# dgelist.general <- DGEList(counts=count.table, group=paste(library.table$LINE, library.table$SEX, sep="_"))
# cat("Normalizing library size by TMM method built into EdgeR.\n\n")
# dgelist.general <- calcNormFactors(dgelist.general, method="TMM")

# cat("Female Library Size Distribution (x1M):\n")
# female.lib.sz <- round(dgelist.general$samples[female.libs,"lib.size"] / (10^6), digits=1)
# cat(min(female.lib.sz), boxplot.stats(female.lib.sz)$stats, max(female.lib.sz), "\n\n")
# cat("Female Library Normalization Factor Distribution:\n")
# female.norm <- round(dgelist.general$samples[female.libs,"norm.factors"], digits=3)
# cat(min(female.norm), boxplot.stats(female.norm)$stats, max(female.norm), "\n\n")

# cat("Male Library Size Distribution (x1M):\n")
# male.lib.sz <- round(dgelist.general$samples[male.libs,"lib.size"] / (10^6), digits=1)
# cat(min(male.lib.sz), boxplot.stats(male.lib.sz)$stats, max(male.lib.sz), "\n\n")
# cat("Male Library Normalization Factor Distribution:\n")
# male.norm <- round(dgelist.general$samples[male.libs,"norm.factors"], digits=3)
# cat(min(male.norm), boxplot.stats(male.norm)$stats, max(male.norm), "\n\n")

# rpm.table <- as.data.frame(cpm(dgelist.general, normalized.lib.sizes=T))

# Drop expression features with var = 0
gene.vars <- apply(rpm.table, 1, var)
keep.genes <- gene.vars != 0
if(any(!keep.genes)) {
  cat("Dropping", sum(!keep.genes), "gene features with 0 variance across both sexes.\n")
  rpm.table <- rpm.table[keep.genes,]
}

# Tried selecting top 10K genes by overall mean expr
# But didn't make much difference and might actually bias away from sex-specific genes
# avg.expr <- apply(rpm.table, 1, mean)
# avg.expr <- sort(avg.expr, decreasing=T)
# top.genes <- names(avg.expr)[1:10000]
# rpm.table <- rpm.table[top.genes,]

# Do log2 transform and transpose
# TO DO: Parameterize the pseudo.count
lib.expr <- t(log2(rpm.table + pseudo.count))


# --- Principle Component Analysis --- #

# Now do PCA on the libraries:
cat("Running PCA on normalized library expression...\n")
lib.pca <- prcomp(lib.expr)
prop.var <- (lib.pca$sdev^2) / sum(lib.pca$sdev^2)

cat("PC1 explains", round(prop.var[1]*100), "% of total variance.\n\n")

female.PC1 <- lib.pca$x[female.libs,1]
male.PC1 <- lib.pca$x[male.libs,1]

cat("Distribution of PC1 values for Female values:\n")
cat(min(female.PC1), boxplot.stats(female.PC1)$stats, max(female.PC1), "\n\n")

cat("Distribution of PC1 values for Male values:\n")
cat(min(male.PC1), boxplot.stats(male.PC1)$stats, max(male.PC1), "\n\n")

# if((min(male.PC1) < max(female.PC1)) & (max(male.PC1) > min(female.PC1))) {
#  cat("WARNING: Male and Female samples overlap on PC1!\n\n")
# }

# Run LDA on PC1 to predict sex:
# (Converting SEX column to boolean, with Male=1)
lda.table <- as.data.frame(cbind(SEX=(library.table$SEX=="M"), PC1=lib.pca$x[,1]))
sex.lda <- lda(SEX ~ PC1, lda.table)
sex.ldap <- predict(sex.lda)
lda.table$SEX.PRED <- sex.ldap$class

# Add SEX.PRED column to library table and convert back to F/M format
library.table[lda.table$SEX.PRED==0,"SEX.PRED"] <- "F"
library.table[lda.table$SEX.PRED==1,"SEX.PRED"] <- "M"

# Also add the probability of the more likely sex
library.table[lda.table$SEX.PRED==0,"PRED.CONF"] <- sex.ldap$posterior[lda.table$SEX.PRED==0,"0"]
library.table[lda.table$SEX.PRED==1,"PRED.CONF"] <- sex.ldap$posterior[lda.table$SEX.PRED==1,"1"]

# Check if any libraries have prediction != label
bad.libs <- union(row.names(library.table)[library.table$SEX != library.table$SEX.PRED], row.names(library.table)[library.table$PRED.CONF < 0.95])
# Also flag any 
if(length(bad.libs) > 0) {
  cat("\n")
  cat(length(bad.libs), "libraries appear to have wrong or low-confidence SEX label:\n\n")
  print(library.table[bad.libs,c("BATCH","SEX","SEX.PRED","PRED.CONF")])
  cat("\n")
}

# Draw first 2 PCs
lib.colors <- rep("black", times=nrow(library.table))
names(lib.colors) <- row.names(library.table)
lib.colors[female.libs] <- "red3"
lib.colors[male.libs] <- "blue3"

pdf(paste0(myParams$qcpath, "libsex_PCA_plot.pdf"))
  plot(lib.pca$x[,1], lib.pca$x[,2], type="p", col=lib.colors, pch=20, xlab="PC1", ylab="PC2", main="PCA of Libraries by Known Gene RPM")
  # Re-draw the misclassified libraries on this plot with different symbol so they stand out
  if(length(bad.libs) > 0) {
    points(lib.pca$x[bad.libs,1], lib.pca$x[bad.libs,2], type="p", col=lib.colors[bad.libs], pch=19, cex=1)
    points(lib.pca$x[bad.libs,1], lib.pca$x[bad.libs,2], type="p", col=lib.colors[bad.libs], pch=1, cex=2)
  }
dev.off()
cat("\n")


# The rest of the code is experimental and should only be explored in interactive mode, so for now the script quits here
cat("Script completed successfully!\n")
print(proc.time())
quit(save="no")


# TO DO: Re-think the remaining part of this
#  - why might it be invalid for RNA-seq, is it just the bigger dynamic range?
# Could also try limiting to the most sex-specific genes in each direction?

# --- Scaled Expr Quantile Analysis --- #

cat("Splitting expression by sex for scaled quantile analysis.\n")
# Separate male and female expression, then scale
female.expr <- lib.expr[female.libs,]
male.expr <- lib.expr[male.libs,]

female.vars <- apply(female.expr, 2, var)
male.vars <- apply(male.expr, 2, var)

if(any(female.vars == 0)) {
  female.keep <- female.vars > 0
  cat("Dropping", sum(!female.keep), "gene features from female expression due to 0 variance.\n")
  female.expr <- female.expr[,female.keep]
}

if(any(male.vars == 0)) {
  male.keep <- male.vars > 0
  cat("Dropping", sum(!male.keep), "gene features from male expression due to 0 variance.\n")
  male.expr <- male.expr[,male.keep]
}

female.expr.scaled <- scale(female.expr)
male.expr.scaled <- scale(male.expr)

stopifnot(sum(is.na(female.expr.scaled))==0)
stopifnot(sum(is.na(male.expr.scaled))==0)

# TO DO: In baseline data, libraries are already in batch order, 
# but might be useful to guarantee batch order and put some marker lines separating batches

# The quantiles to plot:
plot.quant <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)

# Compute what scaled expr value each quantile maps to in each library:
female.quant.table <- apply(female.expr.scaled, 1, quantile, probs=plot.quant)
male.quant.table <- apply(male.expr.scaled, 1, quantile, probs=plot.quant)

# Use boxplot stats to find outliers for the upper-most quantile
female.upper.stats <- boxplot.stats(female.quant.table[nrow(female.quant.table),], coef=3)
female.upper.outliers <- female.upper.stats$out
female.upper.out.x <- which(colnames(female.quant.table) %in% names(female.upper.outliers))
cat("Female outlier libraries based on", row.names(female.quant.table)[nrow(female.quant.table)], "quantile value:\n")
print(female.upper.outliers)
cat("\n")

male.upper.stats <- boxplot.stats(male.quant.table[nrow(male.quant.table),], coef=3)
male.upper.outliers <- male.upper.stats$out
male.upper.out.x <- which(colnames(male.quant.table) %in% names(male.upper.outliers))
cat("Male outlier libraries based on", row.names(male.quant.table)[nrow(male.quant.table)], "quantile value:\n")
print(male.upper.outliers)
cat("\n")

# Line colors:
line.colors <- colorRampPalette(colors=c("blue3","black","red3"))(length(plot.quant))

# Plot Female results
# X positions for each library
lib.x <- 1:ncol(female.quant.table)
pdf(paste0(my.args["QCPATH"], "libsex_Female_scaled_quantiles.pdf"))
matplot(x=lib.x, y=t(female.quant.table), type="l", col=line.colors, lwd=2, lty=1, ylim=c(-6,6), ylab="Scaled Log2 RPM", xlab="Library Number", main="Female Libraries")
abline(h=4, lty=2)
abline(h=-4, lty=2)
points(female.upper.out.x, female.upper.outliers, pch=1, cex=2)
legend("top", row.names(female.quant.table), title="Quantiles", horiz=T, lwd=2, col=line.colors, bty="n")
dev.off()

# Plot Male results
lib.x <- 1:ncol(male.quant.table)
pdf(paste0(my.args["QCPATH"], "libsex_Male_scaled_quantiles.pdf"))
matplot(x=lib.x, y=t(male.quant.table), type="l", col=line.colors, lwd=2, lty=1, ylim=c(-6,6), ylab="Scaled Log2 RPM", xlab="Library Number", main="Male Libraries")
abline(h=4, lty=2)
abline(h=-4, lty=2)
points(male.upper.out.x, male.upper.outliers, pch=1, cex=2)
legend("top", row.names(male.quant.table), title="Quantiles", horiz=T, lwd=2, col=line.colors, bty="n")
dev.off()


cat("Script completed successfully!\n\n")
print(proc.time())
quit("no")

# --- INTERACTIVE ONLY --- #

# Draw several sample vs sample expr plots - these are tailored to this data set, not generalizable!
# Get the top 50 most female and male-specific genes
sex.diff <- apply(lib.expr[female.libs,], 2, mean) - apply(lib.expr[male.libs,], 2, mean)
female.markers <- names(sex.diff)[order(sex.diff, decreasing=T)[1:50]]
male.markers <- names(sex.diff)[order(sex.diff, decreasing=F)[1:50]]
colors <- rep("grey", times=ncol(lib.expr))
names(colors) <- colnames(lib.expr)
colors[female.markers] <- "red"
colors[male.markers] <- "blue"
plot.order <- c(setdiff(colnames(lib.expr), c(female.markers, male.markers)), female.markers, male.markers)
colors <- colors[plot.order]

# This is a negative control where the "bad" male sample is not actually bad
bad.male <- "843_M1_CGTACG_L002_C5UAVANXX"
good.male <- "843_M2_CCAACA_L007_R1_C7G3HANXX"
good.female <- "843_F1_CAGGCG_L006_C6VFHANXX"
male.quant.table[,c(bad.male,good.male)]
# These two samples are both OK, this is what an M vs F should look like
# Shows clear separation of male and female marker genes
plot(lib.expr[good.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# This is the bad sample (X axis) vs replicate of same sex
# In negative control, all sex-markers are correlated, all male sex markers are higher in both samples
plot(lib.expr[bad.male,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)
# There does appear to be some female contamination in this sample based on red dots
# And here is same bad sample (X axis) vs female sample
# Again, in negative control, clear separation of markers
plot(lib.expr[bad.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# Some of the female-specific genes are more highly expr in female but are present in this male sample

# This is a positive example - clear female contamination in 646_M1 lib
bad.male <- "646_M1_GTTTCG_L007_C6VLNANXX"
good.male <- "646_M2_CAAAAG_L002_R1_C7PG0ANXX"
good.female <- "646_F1_CAGGCG_L007_C6VLNANXX"
male.quant.table[,c(bad.male,good.male)]
# These two samples are both OK, this is what an M vs F should look like
# Clear separation of male and female markers
plot(lib.expr[good.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# This is the bad sample (X axis) vs replicate of same sex
# Here we see some high expression of female marker genes
plot(lib.expr[bad.male,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)
# There does appear to be some female contamination in this sample based on red dots
# And here is same bad sample (X axis) vs female sample
# Separation is not as good, many female genes showing up in the main line of points, although shifted up
plot(lib.expr[bad.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# Some of the female-specific genes are more highly expr in female but are present in this male sample

# Example is a bit more complex - looks like there might be mild contamination in both male samples?
# But could also be biological patterns we don't understand yet, I think leave this one be for now
bad.male <- "461_M2_ACTTGA_L001_R1_C8101ANXX"
good.male <- "461_M1_GATCAG_L003_R1_C7G3HANXX"
good.female <- "461_F2_ATCACG_L007_R1_C81A7ANXX"
male.quant.table[,c(bad.male,good.male)]
# Seperation is not as good to begin with
plot(lib.expr[good.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# This is the bad sample (X axis) vs replicate of same sex
# Female marker genes get a bit high, but correlated in both replicates
# This could just be some other biological pattern
plot(lib.expr[bad.male,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)
# There does appear to be some female contamination in this sample based on red dots
# And here is same bad sample (X axis) vs female sample
plot(lib.expr[bad.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# Some of the female-specific genes are more highly expr in female but are present in this male sample

# This was next worst in terms of quantile, but looks fine
bad.male <- "391_M1_GGTAGC_L008_C6VLNANXX"
good.male <- "391_M2_CATGGC_L001_C6VLNANXX"
good.female <- "391_F2_CTCAGA_L007_R1_C81A7ANXX"
male.quant.table[,c(bad.male,good.male)]
# These two samples are both OK, this is what an M vs F should look like
# Good: Clear separation
plot(lib.expr[good.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# This is the bad sample (X axis) vs replicate of same sex
# Looks good
plot(lib.expr[bad.male,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)
# There does appear to be some female contamination in this sample based on red dots
# And here is same bad sample (X axis) vs female sample
# Looks good
plot(lib.expr[bad.male,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)

# Now check female samples - this one looks OK
bad.female <- "801_F1_CATTTT_L002_C5UAVANXX"
good.female <- "801_F2_TACAGC_L004_R1_C81A7ANXX"
good.male <- "801_M1_CAAAAG_L006_C5UAVANXX"
female.quant.table[,c(bad.female,good.female)]
# These two samples are both OK, this is what an M vs F should look like
# Good: Clear separation
plot(lib.expr[good.female,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)
# This is the bad sample (X axis) vs replicate of same sex
# Looks good
plot(lib.expr[bad.female,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# There does appear to be some female contamination in this sample based on red dots
# And here is same bad sample (X axis) vs female sample
# Looks good
plot(lib.expr[bad.female,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)

# This one also looks fine
library.table[library.table$LINE==634,]
bad.female <- "634_F2_CTATAC_L002_C6VLNANXX"
good.female <- "634_F1_TCGGCA_L003_C6VLNANXX"
good.male <- "634_M2_CATGGC_L005_R1_C81A7ANXX"
female.quant.table[,c(bad.female,good.female)]
# These two samples are both OK, this is what an M vs F should look like
# Good: Clear separation
plot(lib.expr[good.female,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)
# This is the bad sample (X axis) vs replicate of same sex
# Looks good
plot(lib.expr[bad.female,plot.order], lib.expr[good.female,plot.order], type="p", pch=20, col=colors)
# There does appear to be some female contamination in this sample based on red dots
# And here is same bad sample (X axis) vs female sample
# Looks good
plot(lib.expr[bad.female,plot.order], lib.expr[good.male,plot.order], type="p", pch=20, col=colors)






