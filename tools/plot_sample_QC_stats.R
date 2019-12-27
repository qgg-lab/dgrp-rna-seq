#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 1/28/15
# Main goal is to visualize the data put together by collect_sample_QC_stats.R and stored in sample_QC_stats.txt
#
# USAGE:
# Rscript plot_sample_QC_stats.R [statfile.txt outputdir]
#
# Default (no additional params): Analyzes QC stats in sample_QC_stats.txt and outputs all figures to sample_QC_figures/
#
# NOTE: This script was originally written for DGRP Baseline analysis, but has now been adapted to be more general
#
# TO DO:
#  - Encapsulate the main code into a "main" function, and only run it if run from command line
#		(In other words, when this script is called from a source() function, it should just load the library of plotting functions
#  - Set the y-limits for each column type automatically, and maybe create command-line params to set manually?
#  - Create violin plots for the main stats (both comparing lanes within batches and across batches)
#  - Add best-fit lines to the technical and biological correlation plots (Maybe move this section to another script?)
#

# setwd("~/Projects/DGRP_Baseline_RNAseq")

library(reshape2)
options(stringsAsFactors=F)

# If two command line parameters are given, change the input table
# and output.dir (But previous settings are the default)
# This will allow updating of the premicrobe plots to upgraded formats for comparison
usageStr <- "USAGE:\nRscript plot_sample_QC_stats.R [statfile.txt outputdir]"
myArgs <- commandArgs(trailingOnly=T)
if(length(myArgs) == 0) {
  # Defaults:
  sample.stat.file <- "sample_QC_stats.txt"
  output.dir <- "sample_QC_figures"
} else if(length(myArgs) == 2) {
  sample.stat.file <- myArgs[1]
  output.dir <- myArgs[2]
} else {
  stop("The script requires exactly 2 parameters, or none for defaults. ", usageStr)
}

# Load the data from stat file
if(file.exists(sample.stat.file)) {
  cat("Loading stats from",sample.stat.file,"\n")
  sample.stat.table <- read.table(sample.stat.file, header=T, sep="\t", quote="")
} else {
  stop("Cannot find ", sample.stat.file, ", either specify correct statfile by command line or run collect_sample_QC_stats.R. ", usageStr)
}

# Compute the fraction that failed to align
# NOTE: "too short" is actually the percentage of reads where the seed region was too short
# The remaining fraction are reads with no seed matching at all, or perhaps other reasons for failure to align
# So "NO.ALIGN" will be the remaining fraction
# "UNALIGNED" will be = "TOO.SHORT" + "NO.ALIGN"
sample.stat.table[,"FRAC.NO.ALIGN"] <- 1 - (sample.stat.table$FRAC.UNIQUE.ALIGN + sample.stat.table$FRAC.MULTI.ALIGN + sample.stat.table$FRAC.TOO.SHORT)
sample.stat.table[,"FRAC.UNALIGNED"] <- sample.stat.table$FRAC.TOO.SHORT + sample.stat.table$FRAC.NO.ALIGN

# All figures will go in this directory
if(!file.exists(output.dir)) {
  dir.create(output.dir, recursive=T)
}

# --- Other Parameters --- #

lane.color.scheme <- c("red1","cyan1","yellow1","purple1","orange1","blue1","green1","violet")
sex.color.scheme <- c(F="red2", M="blue2")
# Number of batches and names will vary by project
# TO DO: Instead of fixed set of colors to use, automatically create a "rainbow" of the necessary number of colors
# TO DO: Create parameters or an additional script to explicitly define the batch structure and corresponding colors
batch.color.scheme <- c("orange1", "forestgreen", "purple", "cyan1", "springgreen2", "blue1", "red1", "slateblue", "gold1", "salmon1", "seagreen1","orangered1")

# --- Plotting Functions --- #

# Core barplot function that:
# A) Optionally sorts the samples by the plotting value column
# B) Optionally colors the bars by sex, lane, or any other parameter
# C) Optionally groups the samples by another column (lane, sex, etc.), and orders the groups
# TO DO: It would be nice if this determined ylim based on the full table, and keeps it constant across multiple 
stats.barplot <- function(
  plot.col,
  scale.values=1,
  order.by.value=F,
  value.decreasing=F,
  group.by.col=NA,
  group.decreasing=F,
  name.col="SAMPLE.NAME",
  stat.table=sample.stat.table,
  color.by.col=NA,  # Column to color by, if NA all bars get same color
  color.scheme="grey70",  # If color.by.col, this should be a vector that maps all possible values in color.by.col to a different color
  split.by.col=NA, # Column to split on for separate plots (each pdf file will have the split column value at end of file name )
  plot.pdf=T,
  pdf.file="",  # If left blank, a name will be automatically generated - this should NOT include the .pdf extension so that the split.by.col value can be automatically added to the end
  pdf.width=9,
  pdf.height=6,
  ylab=plot.col,
  lwd=0.5,
  space=0,
  ...
) {
  # Construct a list of tables to plot - this will contain a single table if split.by.col=NA or if that column has a single value
  # But if not, it will generate multiple pdfs
  stat.table.list <- list()
  # Split by column?
  if(!is.na(split.by.col)) {
    for(split.value in unique(stat.table[,split.by.col])) {
      stat.table.list[[as.character(split.value)]] <- stat.table[stat.table[,split.by.col]==split.value,]
    }
  } else {
    stat.table.list[[1]] <- stat.table
  }
  # Now loop over the sub tables
  for(stat.subtable in stat.table.list) {
    # Order by plot.col?
    if(order.by.value) {
      stat.subtable <- stat.subtable[order(stat.subtable[,plot.col], decreasing=value.decreasing),]
    }
    # Order by group?
    if(!is.na(group.by.col)) {
      stat.subtable <- stat.subtable[order(stat.subtable[,group.by.col], decreasing=group.decreasing),]
    }
    # Now set up bar heights, colors, labels
    bar.heights <- stat.subtable[,plot.col] * scale.values
    if(is.na(color.by.col)) {
      bar.colors <- rep(color.scheme,times=length(bar.heights))
    } else {
      # TO DO: It would be nice if color.scheme could be automatically generated for an arbitrary number of groups
      bar.colors <- color.scheme[stat.subtable[,color.by.col]]
    }
    if(!is.na(name.col)) {
      bar.names <- stat.subtable[,name.col]
    }
    if(plot.pdf) {
      if(pdf.file == "") {
        pdf.file <- paste(output.dir, "/", plot.col, sep="")
        if(!is.na(group.by.col)) {
          pdf.file <- paste(pdf.file, "_by_", group.by.col, sep="")
        }
      }
      sub.pdf.file <- pdf.file
      if(length(stat.table.list) > 1) {
        stopifnot(length(unique(stat.subtable[,split.by.col])) == 1)
        sub.pdf.file <- paste(sub.pdf.file, "_", unique(stat.subtable[,split.by.col]), sep="")
      }
      sub.pdf.file <- paste(sub.pdf.file, ".pdf", sep="")
      # Formatting for presentation slides
      # TO DO: Make these part of the main function paramters
      pdf(sub.pdf.file, width=pdf.width, height=pdf.height, pointsize=20)
      par(mar=c(2.4,3,0.8,0.5), mgp=c(1.6,0.5,0))
    }
    bar.pos <- barplot(bar.heights, col=bar.colors, ylab=ylab, lwd=lwd, space=space, ...)
    if(!is.na(name.col)) {
      mtext(bar.names, side=1, at=bar.pos, las=3, cex=0.5)
    }
    if(plot.pdf) {
      dev.off()
    }
  }
}


# Barplot of a value by lane number
barplot.by.lane <- function(
  plot.col,
  order.by.value=T,
  group.by.col="LANE",
  color.by.col="LANE",
  color.scheme=lane.color.scheme,
  split.by.col="BATCH",
  ...
) {
  stats.barplot(plot.col=plot.col, order.by.value=order.by.value, group.by.col=group.by.col, color.by.col=color.by.col, color.scheme=color.scheme, split.by.col=split.by.col, ...)
}

# Barplot of a value by sex
barplot.by.sex <- function(
  plot.col,
  order.by.value=T,
  group.by.col="SEX",
  color.by.col="SEX",  # Column to color by, if NA all bars get same color
  color.scheme=sex.color.scheme,
  border=NA,
  name.col=NA,
  ...
) {
  stats.barplot(plot.col=plot.col, order.by.value=order.by.value, group.by.col=group.by.col, color.by.col=color.by.col, color.scheme=color.scheme, border=border, name.col=name.col, ...)
}

# Barplot by batch
barplot.by.batch <- function(
  plot.col,
  order.by.value=T,
  group.by.col="BATCH",
  color.by.col="BATCH",
  color.scheme=batch.color.scheme,
  border=NA,
  name.col=NA,
  ...
) {
  stats.barplot(plot.col=plot.col, order.by.value=order.by.value, group.by.col=group.by.col, color.by.col=color.by.col, color.scheme=color.scheme, border=border, name.col=name.col, ...)
}

# Draw the bar plots
barplot.by.lane("FRAC.BASES.TRIMMED", scale.values=100, ylab="% of Bases Trimmed by Cutadapt", ylim=c(0,20))
barplot.by.lane("READ.LEN.MEAN", ylab="Average Read Length", ylim=c(0,125))
barplot.by.lane("READ.LEN.MEDIAN", ylab="Median Read Length", ylim=c(0,125))
barplot.by.lane("FRAC.READS.RIBO", scale.values=100, ylab="% Ribosomal Reads", ylim=c(0,20))
if("FRAC.READS.MICROBE" %in% colnames(sample.stat.table)) {
  barplot.by.lane("FRAC.READS.MICROBE", scale.values=100, ylab="% Microbial Reads", ylim=c(0,32))
}
if("FRAC.READS.REPEAT" %in% colnames(sample.stat.table)) {
  barplot.by.lane("FRAC.READS.REPEAT", scale.values=100, ylab="% Repeat Element Reads", ylim=c(0,10))
}
barplot.by.lane("FILTERED.READS", scale.values=10**-6, ylab="Filtered Read Count (Millions)", ylim=c(0,35))
barplot.by.lane("UNIQUE.ALIGN.READS", scale.values=10**-6, ylab="Uniquely Aligned Reads (Millions)", ylim=c(0,30))
barplot.by.lane("FRAC.UNIQUE.ALIGN", scale.values=100, ylab="% of Uniquely Aligned Reads", ylim=c(0,90))
barplot.by.lane("FRAC.MULTI.ALIGN", scale.values=100, ylab="% of Multiple Aligned Reads", ylim=c(0,30))
barplot.by.lane("FRAC.TOO.SHORT", scale.values=100, ylab="% of Reads Failed Align (Seed Too Short)", ylim=c(0,35))
barplot.by.lane("FRAC.NO.ALIGN", scale.values=100, ylab="% of Reads Failed Align (No Seed/Other)", ylim=c(0,10))
barplot.by.lane("FRAC.UNALIGNED", scale.values=100, ylab="% of Unaligned Reads", ylim=c(0,35))
barplot.by.lane("FRAC.KNOWN.GENES", scale.values=100, ylab="% of Reads from Known Genes", ylim=c(0,100))

# Plot QC stats grouped by SEX, but ONLY if that information is in the QC table in the first place
if("SEX" %in% colnames(sample.stat.table)) {
  barplot.by.sex("FRAC.BASES.TRIMMED", scale.values=100, ylab="% of Bases Trimmed by Cutadapt")
  barplot.by.sex("READ.LEN.MEAN", ylab="Average Read Length")
  barplot.by.sex("READ.LEN.MEDIAN", ylab="Median Read Length")
  barplot.by.sex("FRAC.READS.RIBO", scale.values=100, ylab="% Ribosomal Reads")
  if("FRAC.READS.MICROBE" %in% colnames(sample.stat.table)) {
    barplot.by.sex("FRAC.READS.MICROBE", scale.values=100, ylab="% Microbial Reads")
  }
  if("FRAC.READS.REPEAT" %in% colnames(sample.stat.table)) {
    barplot.by.sex("FRAC.READS.REPEAT", scale.values=100, ylab="% Repeat Element Reads")
  }
  barplot.by.sex("FILTERED.READS", scale.values=10**-6, ylab="Filtered Read Count (Millions)")
  barplot.by.sex("UNIQUE.ALIGN.READS", scale.values=10**-6, ylab="Uniquely Aligned Reads (Millions)")
  barplot.by.sex("FRAC.UNIQUE.ALIGN", scale.values=100, ylab="% of Uniquely Aligned Reads")
  barplot.by.sex("FRAC.MULTI.ALIGN", scale.values=100, ylab="% of Multiple Aligned Reads")
  barplot.by.sex("FRAC.TOO.SHORT", scale.values=100, ylab="% of Reads Failed Align (Seed Too Short)")
  barplot.by.sex("FRAC.NO.ALIGN", scale.values=100, ylab="% of Reads Failed Align (No Seed/Other)")
  barplot.by.sex("FRAC.UNALIGNED", scale.values=100, ylab="% of Unaligned Reads")
  barplot.by.sex("FRAC.KNOWN.GENES", scale.values=100, ylab="% of Reads from Known Genes")
}

# Match batch names to names of batch.color.scheme
batch.names <- unique(sample.stat.table[,"BATCH"])
batch.color.scheme <- batch.color.scheme[1:length(batch.names)]
names(batch.color.scheme) <- batch.names
cat("Color scheme for cross-batch plots:\n")
print(batch.color.scheme)
cat("\n")

# Draw figures that pool everything and groups just by batch 
# (ignore lanes, but split batch 1 into batch1A and batch1B for known bead effect, and split Batch 6 into A and B flow cells)
barplot.by.batch("FRAC.BASES.TRIMMED", scale.values=100, ylab="% of Bases Trimmed by Cutadapt", ylim=c(0,20))
barplot.by.batch("READ.LEN.MEAN", ylab="Average Read Length", ylim=c(0,125))
barplot.by.batch("READ.LEN.MEDIAN", ylab="Median Read Length", ylim=c(0,125))
barplot.by.batch("FRAC.READS.RIBO", scale.values=100, ylab="% Ribosomal Reads", ylim=c(0,20))
if("FRAC.READS.MICROBE" %in% colnames(sample.stat.table)) {
  barplot.by.batch("FRAC.READS.MICROBE", scale.values=100, ylab="% Microbial Reads", ylim=c(0,32))
}
if("FRAC.READS.REPEAT" %in% colnames(sample.stat.table)) {
  barplot.by.batch("FRAC.READS.REPEAT", scale.values=100, ylab="% Repeat Element Reads", ylim=c(0,10))
}
barplot.by.batch("FILTERED.READS", scale.values=10**-6, ylab="Filtered Read Count (Millions)", ylim=c(0,35))
barplot.by.batch("UNIQUE.ALIGN.READS", scale.values=10**-6, ylab="Uniquely Aligned Reads (Millions)", ylim=c(0,30))
barplot.by.batch("FRAC.UNIQUE.ALIGN", scale.values=100, ylab="% of Uniquely Aligned Reads", ylim=c(0,90))
barplot.by.batch("FRAC.MULTI.ALIGN", scale.values=100, ylab="% of Multiple Aligned Reads", ylim=c(0,30))
barplot.by.batch("FRAC.TOO.SHORT", scale.values=100, ylab="% of Reads Failed Align (Seed Too Short)", ylim=c(0,35))
barplot.by.batch("FRAC.NO.ALIGN", scale.values=100, ylab="% of Reads Failed Align (No Seed/Other)", ylim=c(0,10))
barplot.by.batch("FRAC.UNALIGNED", scale.values=100, ylab="% of Unaligned Reads", ylim=c(0,35))
barplot.by.batch("FRAC.KNOWN.GENES", scale.values=100, ylab="% of Reads from Known Genes", ylim=c(0,100))


# -- Draw a figure that piles the counts for same samples on top of each other -- #

# REMOVING THIS PLOT, IT WAS NOT VERY USEFUL
# # In order to draw a barplot where each bar segment is color-coded by the batch it came from
# # need to create a count matrix that is all samples x all batches
# # The cells where that sample was not sequenced in that batch will be 0 and should not show up in the barplot
# 
# # Fix sample names (this could be done at the beginning...)
# sample.stat.table$SAMPLE.NAME <- sub("-","_",sample.stat.table$SAMPLE.NAME,fixed=T)
# 
# # Construct table with UNIQUE.ALIGN.READS for each sample, batch pair
# # This first command should naturally keep samples grouped together by the first batch they were sequenced in
# all.samples <- unique(sample.stat.table$SAMPLE.NAME)
# sample.seq.table <- matrix(data=0, nrow=length(all.samples), ncol=length(batch.color.scheme), dimnames=list(all.samples, names(batch.color.scheme)))
# 
# # Loop over each batch, fill in the values for samples sequenced in that batch
# for(batch in colnames(sample.seq.table)) {
#   batch.orig.rows <- which(sample.stat.table$BATCH == batch)
#   batch.samples <- sample.stat.table$SAMPLE.NAME[batch.orig.rows]
#   sample.seq.table[batch.samples,batch] <- sample.stat.table$UNIQUE.ALIGN.READS[batch.orig.rows]
# }
# 
# Output this table for later
# write.table(cbind(Sample=row.names(sample.seq.table), sample.seq.table), paste0(output.dir, "/sample_batch_unique_aligned_read_counts.txt"), row.names=F, quote=F, sep="\t")
# 
# pdf(paste(output.dir, "/UNIQUE.ALIGN.READS_all_samples_stacked.pdf", sep=""), width=9, height=6)
#   barplot(t(sample.seq.table)*(10**-6), beside=F, space=0, col=batch.color.scheme, ylab="Uniquely Aligned Reads (Millions)", border=NA)
# dev.off()
# 


# -- Draw a figure that sums the counts per sample and draws by-sex plot -- #

# Masking this for now - we don't have both sexes for all projects
# # Construct a new data frame with a single row for each unique sample name, with the sample annotation columns, and the total UNIQUE.ALIGN.READS
# sample.sum.table <- data.frame(SAMPLE.NAME=row.names(sample.seq.table), UNIQUE.ALIGN.READS=apply(sample.seq.table, 1, sum))
# sample.sum.table[,"DGRP.LINE"] <- sub("_[MF][0-9]$","",sample.sum.table$SAMPLE.NAME)
# sample.sum.table[,"SEX"] <- sub("[0-9]$","",sub("^[0-9]+_","",sample.sum.table$SAMPLE.NAME))
# sample.sum.table[,"REP"] <- sub("[0-9]+_[MF]","",sample.sum.table$SAMPLE.NAME)
# 
# # TO DO: Output this table for later?
# 
# barplot.by.sex(stat.table=sample.sum.table, "UNIQUE.ALIGN.READS", scale.values=10**-6, ylab="Uniquely Aligned Reads (Millions)", ylim=c(0,40), pdf.file=paste(output.dir, "/UNIQUE.ALIGN.READS_sample_sum_by_SEX", sep=""), border=NA, name.col=NA)
# 


# --- Replicate Pair Scatter Plots --- #

# This section probably belongs in a separate script, for projects with sufficient numbers of lines
# Masking it for now
# # Build tables of technical replicate pairs
# # Store the results in a list of tables - each table corresponds to a different stat, each table is two columns (row names indicate the pair of samples)
# tech.rep.stats <- c("FRAC.READS.TRIMMED", "FRAC.BASES.TRIMMED", "FRAC.READS.RIBO", "READ.LEN.MEDIAN", "READ.LEN.MEAN", "FRAC.UNIQUE.ALIGN", "FRAC.MULTI.ALIGN", "FRAC.TOO.SHORT", "FRAC.NO.ALIGN", "FRAC.UNALIGNED", "FRAC.KNOWN.GENES")
# if("FRAC.READS.MICROBE" %in% colnames(sample.stat.table)) {
#   tech.rep.stats <- c(tech.rep.stats, "FRAC.READS.MICROBE")
# }
# if("FRAC.READS.REPEAT" %in% colnames(sample.stat.table)) {
#   tech.rep.stats <- c(tech.rep.stats, "FRAC.READS.REPEAT")
# }
# stopifnot(all(tech.rep.stats %in% colnames(sample.stat.table)))
# tech.rep.data <- list()
# for(stat in tech.rep.stats) {
#   tech.rep.data[[stat]] <- data.frame(X=numeric(0), Y=numeric(0))
# }
# next.row <- 1
# # Also keep track of color coding (indicates sex)
# tech.color.scheme <- character(0)
# # Now loop over each unique sample name
# for(sn in unique(sample.stat.table$SAMPLE.NAME)) {
#   # Get rows for all technical replicates
#   sn.rep.table <- sample.stat.table[sample.stat.table$SAMPLE.NAME==sn,]
#   if(nrow(sn.rep.table) > 1) {
#     # Make sure barcode, line, sex, and rep are identical across tech replicates
#     stopifnot(length(unique(sn.rep.table$BARCODE))==1)
#     stopifnot(length(unique(sn.rep.table$DGRP.LINE))==1)
#     stopifnot(length(unique(sn.rep.table$SEX))==1)
#     stopifnot(length(unique(sn.rep.table$REP))==1)
#     # Loop over all replicate pairs
#     for(i in 1:(nrow(sn.rep.table)-1)) {
#       for(j in (i+1):nrow(sn.rep.table)) {
#         ij.name <- paste(sn.rep.table[i,"SAMPLE.ID"], sn.rep.table[j,"SAMPLE.ID"], sep=":")
#         # Assign color based on sex
#         tech.color.scheme[ij.name] <- sex.color.scheme[unique(sn.rep.table$SEX)]
#         # Loop over each stat
#         for(stat in tech.rep.stats) {
#           tech.rep.data[[stat]][next.row,"X"] <- sn.rep.table[i,stat]
#           tech.rep.data[[stat]][next.row,"Y"] <- sn.rep.table[j,stat]
#           row.names(tech.rep.data[[stat]])[next.row] <- ij.name
#         }
#         next.row <- next.row + 1
#       }
#     }
#   }
# }
# stopifnot(length(unique(unlist(lapply(tech.rep.data, nrow))))==1)
# stopifnot(unique(unlist(lapply(tech.rep.data, nrow)))==length(tech.color.scheme))
# 
# # Draw each plot
# for(stat in names(tech.rep.data)) {
#   cat("Plotting",stat,"for all technical replicate pairs.\n")
#   pdf(paste0(output.dir, "/tech_reps_", stat, "_scatter.pdf"))
#     plot(tech.rep.data[[stat]]$X, tech.rep.data[[stat]]$Y, col=tech.color.scheme, main=paste(stat,"of Technical Replicates"), xlab=stat, ylab=stat)
#   dev.off()
#   # Report correlation stats:
#   cat("Pearson Correlation:\n")
#   print(cor.test(tech.rep.data[[stat]]$X, tech.rep.data[[stat]]$Y, method="pearson"))
#   cat("Spearman Correlation:\n")
#   print(cor.test(tech.rep.data[[stat]]$X, tech.rep.data[[stat]]$Y, method="spearman"))
#   cat("\n\n")
# }
# 
# # Check some stats against each other to check for correlations,
# # SPECIFICALLY are the unaligned and too short rates correlated with each other or any other parameter that might give us a clue about what the main cause is?
# for(i in 1:(length(tech.rep.stats)-1)) {
#   for(j in (i+1):length(tech.rep.stats)) {
#     stat.i <- tech.rep.stats[i]
#     stat.j <- tech.rep.stats[j]
#     cat("Plotting",stat.i,"vs",stat.j,"for all individual libraries.\n")
#     pdf(paste0(output.dir, "/", stat.i, "_vs_", stat.j, "_scatter.pdf"))
#       plot(sample.stat.table[,stat.i], sample.stat.table[,stat.j], col=sex.color.scheme[sample.stat.table[,"SEX"]], main="All Individual Libraries", xlab=stat.i, ylab=stat.j)
#     dev.off()
#     cat("Pearson correlation:\n")
#     print(cor.test(sample.stat.table[,stat.i], sample.stat.table[,stat.j], method="pearson"))
#     cat("Spearman Correlation:\n")
#     print(cor.test(sample.stat.table[,stat.i], sample.stat.table[,stat.j], method="spearman"))
#     cat("\n\n")
#   }
# }
# 
# 
# # Collapse the tech replicate table to biological samples and do the same thing for all biological replicate pairs (sum the read counts and recompute the fractions)
# # Technical things like FRAC.READS.RIBO should not be correlated across biological samples
# # But things like FRAC.MULTI.ALIGN and FRAC.UNALIGNED might be
# # NOTE: FRAC.TOO.SHORT counts reads for which the alignable portion was too short (not necessarily the read itself)
# # FRAC.UNALIGNED consists of: reads with too many alignments, reads with too many mismatches, and reads unaligned for "other" reasons
# # May want to break these groups out further...
# sample.stat.table[,"TRIMMED.READS"] <- round(sample.stat.table[,"TOTAL.READS"] * sample.stat.table[,"FRAC.READS.TRIMMED"])
# sample.stat.table[,"RIBO.READS"] <- round(sample.stat.table[,"TOTAL.READS"] * sample.stat.table[,"FRAC.READS.RIBO"])
# if("FRAC.READS.MICROBE" %in% colnames(sample.stat.table)) {
#   sample.stat.table[,"MICROBE.READS"] <- round((sample.stat.table[,"TOTAL.READS"] - (sample.stat.table[,"TRIMMED.READS"] + sample.stat.table[,"RIBO.READS"])) * sample.stat.table[,"FRAC.READS.MICROBE"])
# }
# if("FRAC.READS.REPEAT" %in% colnames(sample.stat.table)) {
#   sample.stat.table[,"REPEAT.READS"] <- round((sample.stat.table[,"TOTAL.READS"] - apply(sample.stat.table[,c("TRIMMED.READS","RIBO.READS","MICROBE.READS")], 1, sum)) * sample.stat.table[,"FRAC.READS.REPEAT"])
# }
# sample.stat.table[,"TOO.SHORT.READS"] <- round(sample.stat.table[,"FILTERED.READS"] * sample.stat.table[,"FRAC.TOO.SHORT"])
# sample.stat.table[,"NO.ALIGN.READS"] <- round(sample.stat.table[,"FILTERED.READS"] * sample.stat.table[,"FRAC.NO.ALIGN"])
# sample.stat.table[,"UNALIGNED.READS"] <- round(sample.stat.table[,"FILTERED.READS"] * sample.stat.table[,"FRAC.UNALIGNED"])
# 
# # Compute sample-level stats using melt - require unique sample name AND barcode just in case
# id.cols <- c("SAMPLE.NAME","BARCODE")
# count.cols <- grep("[.]READS$", colnames(sample.stat.table), value=T)
# lib.qc.melt <- melt(sample.stat.table[,c(id.cols, count.cols)], id.vars=id.cols, measure.vars=count.cols)
# collapsed.qc.data <- dcast(lib.qc.melt, SAMPLE.NAME + BARCODE ~ variable, fun.aggregate=sum)
# # Are all sample names unique?
# if(any(duplicated(collapsed.qc.data$SAMPLE.NAME))) {
#   cat("WARNING: Some sample names occur more than once:", collapsed.qc.data[duplicated(collapsed.qc.data$SAMPLE.NAME),"SAMPLE.NAME"], "\n\n")
# }
# # Parse sample names to get DGRP line, sex, and replicate
# collapsed.qc.data[,"DGRP.LINE"] <- unlist(lapply(strsplit(collapsed.qc.data$SAMPLE.NAME, split="_", fixed=T), "[", 1))
# collapsed.qc.data[,"SEX"] <- as.character(NA)
# collapsed.qc.data[grep("M", collapsed.qc.data$SAMPLE.NAME),"SEX"] <- "M"
# collapsed.qc.data[grep("F", collapsed.qc.data$SAMPLE.NAME),"SEX"] <- "F"
# stopifnot(sum(is.na(collapsed.qc.data$SEX))==0)
# collapsed.qc.data[,"REP"] <- as.integer(unlist(lapply(strsplit(collapsed.qc.data$SAMPLE.NAME, split="[MF]", fixed=F), "[", 2)))
# # Recompute all fractions from summed read counts
# collapsed.qc.data[,"FRAC.READS.TRIMMED"] <- collapsed.qc.data[,"TRIMMED.READS"] / collapsed.qc.data[,"TOTAL.READS"]
# collapsed.qc.data[,"FRAC.READS.RIBO"] <- collapsed.qc.data[,"RIBO.READS"] / collapsed.qc.data[,"TOTAL.READS"]
# if("FRAC.READS.MICROBE" %in% colnames(sample.stat.table)) {
#   collapsed.qc.data[,"FRAC.READS.MICROBE"] <- collapsed.qc.data[,"MICROBE.READS"] / (collapsed.qc.data[,"TOTAL.READS"] - apply(collapsed.qc.data[,c("TRIMMED.READS","RIBO.READS")], 1, sum))
# }
# if("FRAC.READS.REPEAT" %in% colnames(sample.stat.table)) {
#   collapsed.qc.data[,"FRAC.READS.REPEAT"] <- collapsed.qc.data[,"REPEAT.READS"] / (collapsed.qc.data[,"TOTAL.READS"] - apply(collapsed.qc.data[,c("TRIMMED.READS","RIBO.READS","MICROBE.READS")], 1, sum))
# }
# collapsed.qc.data[,"FRAC.UNIQUE.ALIGN"] <- collapsed.qc.data[,"UNIQUE.ALIGN.READS"] / collapsed.qc.data[,"FILTERED.READS"]
# collapsed.qc.data[,"FRAC.MULTI.ALIGN"] <- collapsed.qc.data[,"MULTI.ALIGN.READS"] / collapsed.qc.data[,"FILTERED.READS"]
# collapsed.qc.data[,"FRAC.TOO.SHORT"] <- collapsed.qc.data[,"TOO.SHORT.READS"] / collapsed.qc.data[,"FILTERED.READS"]
# collapsed.qc.data[,"FRAC.NO.ALIGN"] <- collapsed.qc.data[,"NO.ALIGN.READS"] / collapsed.qc.data[,"FILTERED.READS"]
# collapsed.qc.data[,"FRAC.UNALIGNED"] <- collapsed.qc.data[,"UNALIGNED.READS"] / collapsed.qc.data[,"FILTERED.READS"]
# collapsed.qc.data[,"FRAC.KNOWN.GENES"] <- collapsed.qc.data[,"KNOWN.GENE.READS"] / collapsed.qc.data[,"UNIQUE.ALIGN.READS"]
# 
# # Build tables of biological replicate pairs
# # Store the results in a list of tables - each table corresponds to a different stat, each table is two columns (row names indicate the pair of samples)
# bio.rep.stats <- c("FRAC.READS.TRIMMED", "FRAC.READS.RIBO", "FRAC.UNIQUE.ALIGN", "FRAC.MULTI.ALIGN", "FRAC.TOO.SHORT", "FRAC.NO.ALIGN", "FRAC.UNALIGNED", "FRAC.KNOWN.GENES", "TOTAL.READS", "FILTERED.READS", "UNIQUE.ALIGN.READS")
# if("FRAC.READS.MICROBE" %in% colnames(collapsed.qc.data)) {
#   bio.rep.stats <- c(bio.rep.stats, "FRAC.READS.MICROBE")
# }
# if("FRAC.READS.REPEAT" %in% colnames(collapsed.qc.data)) {
#   bio.rep.stats <- c(bio.rep.stats, "FRAC.READS.REPEAT")
# }
# stopifnot(all(bio.rep.stats %in% colnames(collapsed.qc.data)))
# bio.rep.data <- list()
# for(stat in bio.rep.stats) {
#   bio.rep.data[[stat]] <- data.frame(X=numeric(0), Y=numeric(0))
# }
# next.row <- 1
# # Also keep track of color coding (indicates sex)
# bio.color.scheme <- character(0)
# # Now loop over each unique sample name
# for(line in unique(collapsed.qc.data$DGRP.LINE)) {
#   for(sex in c("F","M")) {
#     # Get rows for all technical replicates
#     line.sex.table <- collapsed.qc.data[(collapsed.qc.data$DGRP.LINE==line)&(collapsed.qc.data$SEX==sex),]
#     line.sex.table <- line.sex.table[order(line.sex.table$REP, decreasing=F),]
#     if(nrow(line.sex.table) > 1) {
#       # Loop over all replicate pairs
#       for(i in 1:(nrow(line.sex.table)-1)) {
#         for(j in (i+1):nrow(line.sex.table)) {
#           ij.name <- paste(line.sex.table[i,"SAMPLE.NAME"], line.sex.table[j,"SAMPLE.NAME"], sep=":")
#           # Assign color based on sex
#           bio.color.scheme[ij.name] <- sex.color.scheme[sex]
#           # Loop over each stat
#           for(stat in bio.rep.stats) {
#             bio.rep.data[[stat]][next.row,"X"] <- line.sex.table[i,stat]
#             bio.rep.data[[stat]][next.row,"Y"] <- line.sex.table[j,stat]
#             row.names(bio.rep.data[[stat]])[next.row] <- ij.name
#           }
#           next.row <- next.row + 1
#         }
#       }
#     } else {
#       cat("WARNING: Sample without replication:\n")
#       print(line.sex.table[,c("SAMPLE.NAME","DGRP.LINE","SEX","REP")])
#       cat("\n\n")
#     }
#   }
# }
# stopifnot(length(unique(unlist(lapply(bio.rep.data, nrow))))==1)
# stopifnot(unique(unlist(lapply(bio.rep.data, nrow)))==length(bio.color.scheme))
# cat("\n")
# # TO DO: RESOLVE THE UNPAIRED SAMPLE WARNINGS HERE!
# # I think this is looked at more rigorously in another script though...
# 
# # Draw each plot
# for(stat in names(bio.rep.data)) {
#   cat("Plotting",stat,"for all biological replicate pairs.\n")
#   pdf(paste0(output.dir, "/bio_reps_", stat, "_scatter.pdf"))
#   plot(bio.rep.data[[stat]]$X, bio.rep.data[[stat]]$Y, col=bio.color.scheme, main=paste(stat,"of Biological Replicates"), xlab=stat, ylab=stat)
#   dev.off()
#   # Report correlation stats:
#   cat("Pearson Correlation:\n")
#   print(cor.test(bio.rep.data[[stat]]$X, bio.rep.data[[stat]]$Y, method="pearson"))
#   cat("Spearman Correlation:\n")
#   print(cor.test(bio.rep.data[[stat]]$X, bio.rep.data[[stat]]$Y, method="spearman"))
#   cat("\n\n")
# }
# 

cat("Script completed successfully!\n\n")
print(proc.time())
