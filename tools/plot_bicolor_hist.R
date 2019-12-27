#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 6/13/16
#
# plot_bicolor_hist.R
#
# Goal: Load H^2 or Correlation estimates, draw plot of distributions
# This script was adapted from plot_heritability.R and is designed to replace and extend its functionality
# 
# Usage:
# Rscript plot_bicolor_hist.R [OPTIONS] [DATA=]subdir/data_[H2|traitCorr].txt
#  DATA=    Path to file of data to plot - all H2 or trait correlation columns will be plotted in separate files
#           (If a parameter is given with PARAM= part, it is assumed to be the DATA param)
#  TYPE=    H2 or COR (default is AUTO: inferred from DATA file name)
#  FDR=     Cutoff for determining the "significant" portion of H2 or correlation values, Default=0.05
#  

# Automatically puts all output in the same subdir as the data file, with similar basename
#
# TO DO: Create optional parameters to adjust FDR cutoff.
#     Can automatically detect the input file as the one parameter that does NOT contain an = sign
#


# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

options(stringsAsFactors=F)


# -- Fixed Parameters -- #
# Breaks for the histograms
#   For H2, creates 21 equal sized bins centered on 0, 0.05, ..., 1
#   For Cor, creates 41 equal sized bins centered on -1, -0.95, ..., 0, ..., 1
break.points <- list(
  H2=((0:21)/20)-(1/40),
  COR=((0:41)/20)-(41/40)
)
show.labels <- list(
  H2=c("0","0.25","0.5","0.75","1"),
  COR=c("-1","-0.75","-0.5","-0.25","0","0.25","0.5","0.75","1")
)
# Color scheme - different color schemes for each sex
colors <- list('F'=c(SIG="red3", NS="white"), 'M'=c(SIG="blue3", NS="white"), 'Pooled'=c(SIG="purple3", NS="white"))

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript plot_bicolor_hist.R [OPTIONS] [DATA=]subdir/data_[model_results|traitCorr].txt"
my.args <- commandArgs(trailingOnly=T)
# INTERACTIVE TESTING:
# my.args <- c("transposons_class/combined_transposon_filtered_rpm_WolAdj_model_results.txt")
# my.args <- c("transposons_class/combined_transposon_filtered_rpm_WolAdj_Pooled_model_results.txt")
# my.args <- c("microbiome/combined_microbe_filtered_rpm_model_results.txt")
# my.args <- c("known_novel_genes/combined_samples_split_known_fpkm_Wol_H2.txt")
# my.args <- c("known_novel_genes/traitCorr/combined_samples_known_novel_fpkm_Wol_F_traitCorr.txt")
if(length(my.args) == 0) {
  stop("Requires at least one argument. ", usageStr)
} else if(sum(!grepl("=", my.args, fixed=T)) > 1) {
  stop("Too many unnamed parameters. ", usageStr)
} else {
  unnamed.arg <- grep("=", my.args, invert=T, fixed=T)
  my.args[unnamed.arg] <- paste0("DATA=", my.args[unnamed.arg])
  cat("Running with params:",my.args,"\n")
}

# Parse param names and values
argSplit=strsplit(my.args, split="=", fixed=T)
my.args <- lapply(argSplit, "[", 2)
names(my.args) <- unlist(lapply(argSplit, "[", 1))

# Make sure EXPR param exists
if(!("DATA" %in% names(my.args))) {
  stop("Missing DATA parameter. ", usageStr)
}

# Fill in default value for optional params
setDefault <- function(param, default) {
  if(!(param %in% names(my.args))) {
    cat("Setting ",param,"=",default," by default\n", sep="")
    my.args[[param]] <<- default
  }
}

setDefault("TYPE", "AUTO")

setDefault("FDR", 0.05)
if(!is.numeric(my.args$FDR)) {
  my.args$FDR <- as.numeric(my.args$FDR)
}
if(is.na(my.args$FDR)) {
  cat("WARNING: Invalid FDR, defaulting to 0.05")
  my.args$FDR <- 0.05
}

# Extract path to EXPR input file, use it as the output dir
my.args$OUTDIR <- paste0(dirname(my.args$DATA), "/")

# Base name of COUNTS file comes from shaving off subdirectory part of path, and trailing ".txt"
my.args$BASENAME <- sub("[.]txt$", "", basename(my.args$DATA))

# Determine the type from BASENAME
if(my.args$TYPE=="AUTO") {
  if(grepl("[_]H2$", my.args$BASENAME) | grepl("[_]model[_]results$", my.args$BASENAME)) {
    my.args$TYPE <- "H2"
  } else if(grepl("[_]traitCorr$", my.args$BASENAME)) {
    my.args$TYPE <- "COR"
  } else {
    stop("Could not determine input data type, use TYPE=[H2|COR] to specify correct type.")
  }
}

# Iff TYPE=COR, define a SEX param as well (otherwise this will vary by column)
if(my.args$TYPE=="COR") {
  if(grepl("[_]F[_]", my.args$BASENAME)) {
    my.args$SEX <- "F"
  } else if(grepl("[_]M[_]", my.args$BASENAME)) {
    my.args$SEX <- "M"
  } else {
    stop("Could not determine sex for TYPE=COR, file name should contain \"_F_\" or \"_M_\" to denote this.")
  }
}

cat("\n")


# --- LOAD INPUT --- #

# Load the data table
data.table <- read.table(my.args$DATA, header=T, sep="\t", row.names=1)

# Determine the value and FDR columns
if(my.args$TYPE=="H2") {
  value.cols <- colnames(data.table)[grepl("[.]H2$", colnames(data.table)) | grepl("[.]Perc[.]Var$", colnames(data.table))]
  names(value.cols) <- sub("[.]H2$", "", value.cols)
  names(value.cols) <- sub("[.]Perc[.]Var$", "", names(value.cols))
} else {
  value.cols <- colnames(data.table)[grepl("[.]cor$", colnames(data.table))]
  names(value.cols) <- sub("[.]cor$", "", value.cols)
}

fdr.cols <- grep("[.][Pp][Vv]al[.][BF][^.]*$", colnames(data.table), value=T)
names(fdr.cols) <- sub("[.][Pp][Vv]al[.][BF][^.]*$", "", fdr.cols)
# Drop the ".Line" part of these columns names
# TO DO: For Pooled Models need to combine Line and SexLine p-values here
names(fdr.cols) <- sub("[.]L[iI][nN][eE]$", "", names(fdr.cols))

if(length(fdr.cols) > length(value.cols)) {
  # Special case for Pooled model - need to combine LINE and SEX.LINE p-values
  if("Pooled.SEX.LINE.PVal.FDR" %in% fdr.cols) {
    data.table[,"Pooled.LINE.PVal.FDR"] <- apply(data.table[,c("Pooled.LINE.PVal.FDR","Pooled.SEX.LINE.PVal.FDR")], 1, min)
    fdr.cols <- fdr.cols[fdr.cols != "Pooled.SEX.LINE.PVal.FDR"]
  }
}

stopifnot(length(value.cols) == length(fdr.cols))
stopifnot(all(names(value.cols) == names(fdr.cols)))

cat("Plotting", length(value.cols), my.args$TYPE, "columns each containing", nrow(data.table), "features.\n\n")

# Loop over each sample subset - generate histogram bar heights
hist.bars <- list()
for(j in names(value.cols)) {
  # Pull any extreme values into range
  below.range.vals <- data.table[,value.cols[j]] < min(break.points[[my.args$TYPE]])
  above.range.vals <- data.table[,value.cols[j]] > max(break.points[[my.args$TYPE]])
  if(any(below.range.vals)) {
    cat("WARNING: The following", value.cols[j], "values are below standard histogram range and will be counted in lowest bar:\n")
    print(cbind(row.names(data.table)[below.range.vals], data.table[below.range.vals,value.cols[j]]))
    cat("\n")
    data.table[below.range.vals,value.cols[j]] <- min(break.points[[my.args$TYPE]])
  }
  if(any(above.range.vals)) {
    cat("WARNING: The following", value.cols[j], "values are above standard histograms range and will be counted in highest bar:\n")
    print(cbind(row.names(data.table)[above.range.vals], data.table[above.range.vals,value.cols[j]]))
    cat("\n")
    data.table[above.range.vals,value.cols[j]] <- max(break.points[[my.args$TYPE]])
  }
  
  # Split out the H2 values for significant features and everything else
  j.signif <- (data.table[,fdr.cols[j]] <= my.args$FDR)
  j.signif.values <- data.table[j.signif, value.cols[j]]
  j.bg.values <- data.table[!j.signif, value.cols[j]]
  j.signif.hist <- hist(j.signif.values, breaks=break.points[[my.args$TYPE]], plot=F)
  j.bg.hist <- hist(j.bg.values, breaks=break.points[[my.args$TYPE]], plot=F)
  stopifnot(j.signif.hist$mids == j.bg.hist$mids)
  j.hist <- data.frame(row.names=j.signif.hist$mids, 'Signif'=j.signif.hist$counts, 'Not Signif'=j.bg.hist$counts)
  hist.bars[[j]] <- t(j.hist)
}

# Determine y.max for each histogram
y.max <- round(unlist(lapply(hist.bars, function(x){max(apply(x, 2, sum))}))*1.05)
# For TYPE=="H2", match the Wol and regular H2 column max values
if(my.args$TYPE=="H2") {
  if(all(c("F.Wol","M.Wol") %in% names(y.max))) {
    y.max[c("F.Wol","M.Wol")] <- max(y.max[c("F.Wol","M.Wol")])
  }
  y.max[c("F","M")] <- max(y.max[c("F","M")])
}

# TO DO: For COR TYPE, should we make all y.max values match up?

# Loop over each sample subset - draw plots
for(j in names(value.cols)) {
  # Determine sex for this column
  if(my.args$TYPE=="COR") {
    j.sex <- my.args$SEX
  } else if(grepl("^F", j)) {
    j.sex <- "F"
  } else if(grepl("^M", j)) {
    j.sex <- "M"
  } else if(grepl("^Pooled", j)) {
    j.sex <- "Pooled"
  } else {
    stop("Could not determine sex for column ", value.cols[j])
  }
  
  # Determine appropriate xlab and main
  if(my.args$TYPE=="COR") {
    j.xlab <- "Spearman Correlation"
    j.main <- j
  } else if(grepl("[.]Wol$", j)) {
    j.xlab <- "% Var from Wolbachia"
    j.main <- "Wolbachia Effects"
  } else {
    j.xlab <- expression(paste('Estimated ', italic(H^2), sep=''))
    j.main <- "Heritability"
  }
  if(j.sex=="F") {
    j.main <- paste("Female",j.main)
  } else if(j.sex=="M") {
    j.main <- paste("Male",j.main)
  } else if(j.sex=="Pooled") {
    j.main <- paste("Pooled Sex",j.main)
  }
  
  # Determine color scheme based on sex
  j.colors <- colors[[j.sex]]
  
  # Match up labels to the columns they go with
  # (This is b/c of dumb precision errors)
  j.labels <- show.labels[[my.args$TYPE]]
  j.label.cols <- which(round(as.numeric(colnames(hist.bars[[j]])), digits=2) %in% round(as.numeric(show.labels[[my.args$TYPE]]), digits=2))
  
  pdf.file <- paste0(my.args["OUTDIR"], my.args["BASENAME"], "_", gsub("[.]", "", j), "_hist.pdf")
  cat("Drawing", pdf.file, "\n")
  pdf(pdf.file, width=5.5, height=3.5, pointsize=16)
  par(mgp=c(1.4,0.3,0), mar=c(2.5,3,1,0.2))
  bar.centers <- barplot(hist.bars[[j]], col=j.colors, ylim=c(0,y.max[j]), xlab=j.xlab, ylab="Number of Features", main=j.main, tcl=-0.25, xaxt="n")
  
  axis(side=1, at=bar.centers[j.label.cols], labels=j.labels, tick=F)
  # Getting rid of this legend b/c it just gets in the way
  # legend("topright", paste0(my.args$FDR*100, "% FDR"), fill=j.colors["SIG"])
  dev.off()
}


cat("\n\nScript completed successfully!\n\n")
print(proc.time())
quit("no")
