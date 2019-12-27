#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 12/31/16
#
# Quick script to plot the results from genotype_samples.R
# Output goes in sample_QC_figures
#
# USAGE: genotype_validation_plots.R BATCH1 [BATCH2 ...]
#

#
# TO DO LIST:
# Output of samples with LER <= 1 should show the relevant columns
# Samples with LER = 0 (or NA) should be excluded for plotting
# (These are cases where there were no differentiating SNPs)
#

options(stringsAsFactors=F)

# Take the list of batches to analyze from command line params
usageStr <- "USAGE:\nRscript genotype_validation_plots.R BATCH1 [BATCH2 ...]"
myArgs <- commandArgs(trailingOnly=T)
# Interactive code testing: 
# myArgs <- c("batch_141031","batch_150318","batch_150427","batch_150527","batch_150707","batch_150917")
batches <- myArgs

if(length(batches) == 0) {
  # Can't run without at least one batch!
  stop("Must specify one or more batches to verify genotypes. ", usageStr)
}

# Remove trailing .txt from these files (just in case)
batches <- sub("[.]txt$", "", batches)

cat("Analyzing",length(batches),"batches:", batches, "\n\n")

# Load each batch genotype_check.txt file, combine into single table
check.table <- NULL
for(batch in batches) {
  batch.check.file <- paste0("allele_counts/", batch, "/genotype_check.txt")
  if(file.exists(batch.check.file)) {
    batch.check <- read.table(batch.check.file, header=T, sep="\t")
    cat("Loaded data for", nrow(batch.check), "samples from", batch.check.file, "\n")
    check.table <- rbind(check.table, batch.check)
  } else {
    cat("WARNING: Could not find", batch.check.file, "\n")
  }
}

cat("Loaded genotype validation data for", nrow(check.table), "samples total.\n")
cat(sum(check.table$FLAG=="OK"), "samples have correct genotype label.\n")
cat(sum(check.table$FLAG=="WARN"), "samples have a WARN flag.\n")
cat(sum(check.table$FLAG=="ERR"), "samples have an ERR flag.\n")
cat("\n")

cat("Samples with WARN/ERR flags:\n")
print(check.table[c(which(check.table$FLAG=="WARN"),which(check.table$FLAG=="ERR")),])

# --- PLOT FIGURES --- #
output.dir <- "sample_QC_figures/"
if(!file.exists(output.dir)) {
  dir.create(output.dir, recursive=T)
}

# Figure 1 - Show distribution of Line Mean Error for the annotated line (ANNOT.ERR)
# Show distribution for the good samples (ANNOT.LINE==BEST.LINE), mark the rest with red lines
lme.err.samples <- which(check.table$ANNOT.LINE != check.table$BEST.LINE)
cat(length(lme.err.samples), "samples have ANNOT.LINE != BEST.LINE:\n")
print(check.table[lme.err.samples,c("SAMPLE","FLAG","BEST.LINE","BEST.ERR","ANNOT.LINE","ANNOT.ERR")])
lme.flag.samples <- setdiff(which(check.table$FLAG != "OK"),lme.err.samples)
cat(length(lme.flag.samples), "samples are flagged for other reasons:\n")
print(check.table[lme.flag.samples,c("SAMPLE","FLAG","BEST.LINE","BEST.ERR","ANNOT.LINE","ANNOT.ERR")])
lme.ok.samples <- setdiff(1:nrow(check.table), c(lme.err.samples,lme.flag.samples))
lme.density <- density(check.table[lme.ok.samples,"ANNOT.ERR"])
xlim <- c(0,max(check.table$ANNOT.ERR)*1.05)
ylim <- range(0,max(lme.density$y)*1.05)

# Write out to a PDF in output.dir
pdf(paste0(output.dir, "genotype_annot_lme_boxplot.pdf"), height=2.5, width=5, pointsize=16)
  par(mgp=c(1.8,0.6,0), mar=c(3,0.5,0.5,0.5))
  boxplot(check.table[lme.ok.samples,"ANNOT.ERR"], horizontal=T, xlab="Line Mean Error", ylim=xlim, outpch=20)
  points(x=check.table[lme.flag.samples,"ANNOT.ERR"], y=rep(1,times=length(lme.flag.samples)), pch=4, col="orange3", cex=1.2, lwd=3)  
  points(x=check.table[lme.err.samples,"ANNOT.ERR"], y=rep(1,times=length(lme.err.samples)), pch=4, col="red3", cex=1.2, lwd=3)
dev.off()

# Density version doesn't look good
  
# plot(lme.density$x, lme.density$y, type="l", xlim=xlim, ylim=ylim, xlab="Line Mean Error", ylab="Probability Density", lwd=2)
# for(i in lme.flag.samples) {
#  abline(v=check.table[i,"ANNOT.ERR"], col="orange3", lwd=2)
# }
# for(i in lme.err.samples) {
#  abline(v=check.table[i,"ANNOT.ERR"], col="red3", lwd=2)
# }


# Figure 2 - Show distribution of Line Error Ratio (annotated line vs next best)
# Show distribution for the good samples (LER > 1), mark the rest with red lines
ler.err.samples <- which(check.table$RATIO.ERR <= 1)
cat(length(ler.err.samples), "samples have RATIO.ERR <= 1:\n")
print(check.table[ler.err.samples,c("SAMPLE","FLAG","BEST.LINE","BEST.ERR","ANNOT.LINE","ANNOT.ERR")])
# Not flagging samples that have error flag for other reasons here
ler.flag.samples <- c()
# ler.flag.samples <- setdiff(which(check.table$FLAG != "OK"),ler.err.samples)
# cat(length(ler.flag.samples), "samples are flagged for other reasons:\n")
# print(check.table[ler.flag.samples,c("SAMPLE","FLAG","BEST.LINE","BEST.ERR","ANNOT.LINE","ANNOT.ERR")])
ler.ok.samples <- setdiff(1:nrow(check.table), c(ler.err.samples,ler.flag.samples))
ler.density <- density(check.table[ler.ok.samples,"RATIO.ERR"])
xlim <- c(log10(min(check.table$RATIO.ERR))*1.05,log10(max(check.table$RATIO.ERR))*1.05)
ylim <- range(0,max(ler.density$y)*1.05)

# Write out to a PDF in output.dir
pdf(paste0(output.dir, "genotype_annot_ler_boxplot.pdf"), height=2.5, width=5, pointsize=16)
par(mgp=c(1.8,0.6,0), mar=c(3,0.5,0.5,0.5))
boxplot(check.table[,"RATIO.ERR"], horizontal=T, xlab="Line Error Ratio", outpch=20, log="x", axt="n")
# points(x=check.table[ler.flag.samples,"RATIO.ERR"], y=rep(1,times=length(ler.flag.samples)), pch=4, col="orange3", cex=1.2, lwd=3)  
points(x=check.table[ler.err.samples,"RATIO.ERR"], y=rep(1,times=length(ler.err.samples)), pch=4, col="red3", cex=1.2, lwd=3)
abline(v=1, lty=2)
dev.off()

