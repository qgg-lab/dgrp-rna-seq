#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 1/20/17
#
# filter_line_means.R
#
# The goal of this script is to:
# 1) Filter a line means table based on FLAG columns and average expression (optional)
# 2) Further filter based on H2 and/or FDR (optional - requires providing the right model results table and which gender the line means table is)
# 3) Run hclust on line means after all gene filtering is done, optionally remove outlier lines
# 4) Reformatting for downstream analysis:
#     MMC: replaces commas in row names, then output as CSV
#       TO DO: Should automatically mask meta-features here?
#     eQTL: reformats the same way as plink_setup.R (transposes matrix, adds FID column, optionally permutes)
# This will replace some of these functions (step 4) in the current MMC script and plink scripts
#
# Usage:
# Rscript filter_line_means.R [OPTIONS] EXPR=line_means.txt MODE=[mmc|eqtl]
#  EXPR=  The file containing line mean expression values to filter and prep for downstream analysis
#  MODE=  mmc|eqtl|wgcna (default is mmc)
#  OUTDIR= Directory to output, either specify full path or it will be inferred from path to EXPR file
#         Generally, it takes the path portion of EXPR, drops genVar/, and then adds mmc/ (MODE=mmc), plink/ (MODE=eqtl), or wgcna/ (MODE=wgcna)
#  FILTER=  Flag(s) to filter genes on. Multiple flags should be comma-seperated, no spaces
#           Default is to leave it blank, no filtering on flag column
#           FILTER=RARE filters out all genes with "RARE" in FLAG column of line means table
#  AVGEXPR= Threshold for average expression, in addition to FLAG filtering of line means
#  H2FILE=  The .._model_results.txt file to use for H2 and FDR thresholds
#  SEX=   M/F/Pooled - tells which columns in H2FILE to use.  Tries to infer from line_means.txt file if not specified.
#  FDR=   Max FDR cut-off for which genes/features to include, default = none
#  FDRCOL=  One or more model terms to apply the FDR cutoff on, default = LINE
#           Seperate multiple terms with a comma, no spaces
#  H2=    Minimum H2 cut-off for which genes/features to include, default = none
#  VAR=   Minimum variance cut-off for which genes/features to include, default = none
#  HCUT=  Cut-off point to remove outlier lines based on hclust results, default = none
#  PERM=  eQTL only - after outputting the main table, create this many permuted versions (default=0)
#

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

options(stringsAsFactors=F)

# -- Process command-line parameters -- #
usageStr="USAGE:\nRscript filter_line_means.R [OPTIONS] EXPR=line_means.txt MODE=[mmc|eqtl|wgcna]"
my.args <- commandArgs(trailingOnly=T)
# TEMP TESTING:
# MMC:  my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_F_line_means.txt", "H2FILE=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_model_results.txt", "FILTER=RARE,LOW", "AVGEXPR=0", "H2=0.5", "FDR=0.001")
# MMC Pooled:   my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_Pooled_line_means.txt", "H2FILE=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_Pooled_model_results.txt", "FILTER=RARE,LOW", "AVGEXPR=0", "H2=0.5", "FDR=0.05", "FDRCOL=SEX.LINE") 
# EQTL: my.args <- c("EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_eQTL_F_line_means.txt", "H2FILE=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_eQTL_model_results.txt", "FILTER=RARE,LOW", "MODE=eqtl", "PERM=10")
# WGCNA: my.args <- c("MODE=wgcna","EXPR=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_F_line_means.txt", "H2FILE=known_all_novel_genes/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_model_results.txt", "FILTER=RARE,LOW", "AVGEXPR=0", "H2=0.5", "FDR=0.001")

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

# MODE=  mmc|eqtl|wgcna (default is mmc)
# Default to mmc if not specified, must be mmc or eqtl
setDefault("MODE","mmc")
if(!(my.args$MODE %in% c("mmc","eqtl","wgcna"))) {
  stop(my.args$MODE, " is not a recognized MODE, must be mmc or eqtl. ", usageStr)
}

#  OUTDIR= Directory to output, either specify full path or it will be inferred from path to EXPR file
#         Generally, it takes the path portion of EXPR, drops genVar/, and then adds mmc/ (MODE=mmc) or plink/ (MODE=eqtl)
# 
# Determine the default OUTDIR based on EXPR and MODE...
defOutdir <- paste0(dirname(my.args$EXPR), "/")
defOutdir <- sub("genVar/", "", defOutdir)
if(my.args$MODE=="mmc") {
  defOutdir <- paste0(defOutdir, "mmc/")
} else if(my.args$MODE=="eqtl") {
  defOutdir <- paste0(defOutdir, "plink/")
} else {
  defOutdir <- paste0(defOutdir, "wgcna/")
}
setDefault("OUTDIR",defOutdir)
# Make sure this subdir exists
if(!file.exists(my.args$OUTDIR)) {
  dir.create(path=my.args$OUTDIR, recursive=T, showWarnings = F)
}

#  FILTER=  Flag(s) to filter genes on. Multiple flags should be comma-seperated, no spaces
#           Default is to leave it blank, no filtering on flag column
#           FILTER=RARE filters out all genes with "RARE" in FLAG column of line means table
setDefault("FILTER","")

#  AVGEXPR= Threshold for average expression, in addition to FLAG filtering of line means
setDefault("AVGEXPR",-Inf)
if(class(my.args$AVGEXPR) != "numeric") {
  my.args$AVGEXPR <- as.numeric(my.args$AVGEXPR)
}

#  H2FILE=  The .._model_results.txt file to use for H2 and FDR thresholds
setDefault("H2FILE","")

#  SEX=   M/F/Pooled - tells which columns in H2FILE to use.  Tries to infer from line_means.txt file if not specified.
inferredSex <- regmatches(my.args$EXPR, regexec(".*[_]([A-Za-z]+)[_]line[_][md][ei][af][nf]s[.]txt", my.args$EXPR))[[1]][2]
if((length(inferredSex) == 0) | is.na(inferredSex)) {inferredSex <- ""}
if(inferredSex == "SexAvg") {inferredSex <- "Pooled"}
if(!(inferredSex %in% c("M","F","Pooled"))) {inferredSex <- ""}
# Change SexAvg to Pooled, set to "" if not one of the main values
setDefault("SEX",inferredSex)

# FDR=   Max FDR cut-off for which genes/features to include, default = none
setDefault("FDR",Inf)
if(class(my.args$FDR) != "numeric") {
  my.args$FDR <- as.numeric(my.args$FDR)
}

setDefault("FDRCOL","LINE")

# H2=    Minimum H2 cut-off for which genes/features to include, default = none
setDefault("H2",-Inf)
if(class(my.args$H2) != "numeric") {
  my.args$H2 <- as.numeric(my.args$H2)
}

# If FDR or H2 is > Inf, but no H2FILE specified, quit with error
if(((my.args$FDR < Inf) | (my.args$H2 > -Inf)) & (my.args$H2FILE == "")) {
  stop("Must specify H2FILE in order to set FDR or H2 threshold.")
}

# VAR=  Minimum variance for which genes/features to include, default = none
setDefault("VAR",-Inf)
if(class(my.args$VAR) != "numeric") {
  my.args$VAR <- as.numeric(my.args$VAR)
}

#  HCUT=  Cut-off point to remove outlier lines based on hclust results, default = none
setDefault("HCUT", Inf)

# PERM=  eQTL only - after outputting the main table, create this many permuted versions (default=0)
setDefault("PERM",as.integer(0))
if(class(my.args$PERM) != "integer") {
  my.args$PERM <- as.integer(my.args$PERM)
}

# Set RNG (could parameterize the seed here?)
my.seed <- as.integer(8210678)
cat("Seeding RNG with", my.seed, "\n")
set.seed(my.seed)

cat("\n")


# --- Load input files --- #

cat("Loading expression data from:", my.args$EXPR, "\n")
expr.table <- read.table(my.args$EXPR, header=T, row.names=1, sep="\t", check.names=F)
cat("Loaded expression for", nrow(expr.table), "across", length(setdiff(colnames(expr.table), "FLAG")), "conditions.\n")

# If H2FILE specified, load it now
if(my.args$H2FILE != "") {
  cat("Loading model results from:", my.args$H2FILE, "\n")
  h2.table <- read.table(my.args$H2FILE, header=T, row.names=1, sep="\t")
  cat("Loaded model results for", nrow(h2.table), "features.\n")
  # Make sure there are model results for all the features in expr.table
  if(!all(row.names(expr.table) %in% row.names(h2.table))) {
    stop("Specified H2FILE is missing features in EXPR.")
  }
  h2.table <- h2.table[row.names(expr.table),]
  # Figure out the appropriate FDR column(s)
  fdr.cols <- paste0(unlist(strsplit(my.args$FDRCOL, split=",")), ".PVal.FDR")
  # H2 col is always just H2 for now
  h2.col <- "H2"
  if(my.args$SEX != "") {
    fdr.cols <- paste(my.args$SEX, fdr.cols, sep=".")
    h2.col <- paste(my.args$SEX, h2.col, sep=".")
  }
  # Make sure are both in this table, then reduce to just those columns
  if(!all(fdr.cols %in% colnames(h2.table))) {
    stop("H2FILE is missing ", setdiff(fdr.cols,colnames(h2.table)))
  }
  if(!(h2.col %in% colnames(h2.table))) {
    stop("H2FILE is missing ", h2.col)
  }
  h2.table <- h2.table[,c(h2.col,fdr.cols)]
} else {
  h2.table <- NULL
  h2.col <- NULL
  fdr.cols <- NULL
}


# --- Filter based on user-specified parameters --- #

#  FILTER=  Flag(s) to filter genes on. Multiple flags should be comma-seperated, no spaces
#           Default is to leave it blank, no filtering on flag column
#           FILTER=RARE filters out all genes with "RARE" in FLAG column of line means table
#  AVGEXPR= Threshold for average expression, in addition to FLAG filtering of line means
#  FDR=   Max FDR cut-off for which genes/features to include, default = none
#  H2=    Minimum H2 cut-off for which genes/features to include, default = none
#  HCUT=  Cut point for hclust results to remove outlier lines

# Filter a line means table based on FLAG columns
if(my.args$FILTER != "") {
  if(!("FLAG" %in% colnames(expr.table))) {
    stop("No FLAG column in EXPR table, cannot FILTER out ", my.args$FILTER, " features.")
  }
  cat("Filtering by FLAG column, total feature counts:\n")
  print(table(expr.table$FLAG))
  remove.flags <- unlist(strsplit(my.args$FILTER,split=","))
  cat("Removing features with the following flags:", remove.flags, "\n")
  # Also drop the FLAG column here
  expr.table <- expr.table[!(expr.table$FLAG %in% remove.flags),setdiff(colnames(expr.table),"FLAG")]
  cat(nrow(expr.table), "features remain.\n\n")
} else {
  # Just drop the FLAG column if specified
  if("FLAG" %in% colnames(expr.table)) {
    cat("No specified FLAG filter, removing FLAG column.\n")
    expr.table <- expr.table[,setdiff(colnames(expr.table),"FLAG")]
  }
}

# Filter by average expression (optional)
if(my.args$AVGEXPR > -Inf) {
  cat("Dropping features with average expression <", my.args$AVGEXPR, "\n")
  
  # Compute avg expression
  avg.expr <- apply(expr.table, 1, mean)
  
  # Determine which features are < AVGEXPR cutoff
  above.cutoff <- avg.expr >= my.args$AVGEXPR
  
  cat(sum(!above.cutoff), "features are below threshold.\n")
  expr.table <- expr.table[above.cutoff,]
  cat(nrow(expr.table), "features remain.\n\n")
}

# Filter by H2 cutoff
if(my.args$H2 > -Inf) {
  cat("Dropping features with",h2.col,"<", my.args$H2, "\n")
  
  # Filter h2.table down to the rows still in expr.table
  h2.table <- h2.table[row.names(expr.table),]
  
  above.cutoff <- h2.table[,h2.col] >= my.args$H2
  
  cat(sum(!above.cutoff), "features are below threshold.\n")
  expr.table <- expr.table[above.cutoff,]
  cat(nrow(expr.table), "features remain.\n\n")
}

# Filter by FDR cutoff
if(my.args$FDR < Inf) {
  cat("Dropping features with",paste(fdr.cols, collapse=" & "),">", my.args$FDR, "\n")
  
  # Filter h2.table down to the rows still in expr.table
  h2.table <- h2.table[row.names(expr.table),]
  
  if(length(fdr.cols) > 1) {
    min.fdr <- apply(h2.table[,fdr.cols], 1, min)
  } else {
    min.fdr <- h2.table[,fdr.cols]
  }
  
  below.cutoff <- min.fdr <= my.args$FDR
  
  cat(sum(!below.cutoff), "features are above treshold.\n")
  expr.table <- expr.table[below.cutoff,]
  cat(nrow(expr.table), "features remain.\n\n")
}

# Filter by VAR cutoff
if(my.args$VAR > -Inf) {
  cat("Dropping features with expression variance <", my.args$VAR, "\n")
  
  # Compute expr variance of each feature
  expr.var <- apply(expr.table, 1, var)
  
  above.cutoff <- expr.var >= my.args$VAR
  
  cat(sum(!above.cutoff), "features are below threshold.\n")
  expr.table <- expr.table[above.cutoff,]
  cat(nrow(expr.table), "features remain.\n\n")
}

# Run hclust on line means using remaining features
# This is done regardless of whether HCUT is specified so the user can see if there are any potential problems
lineTree <- hclust(dist(t(expr.table)), method="average")
plot.file <- basename(my.args$EXPR)
plot.file <- sub("[_]line[_][md][ei][af][nf]s[.]txt$", "", plot.file)
plot.file <- paste0(my.args$OUTDIR, plot.file, "_line_hclust.pdf")
pdf(plot.file, width=12, height=9)
par(cex=0.6, mar=c(0,4,2,0))
plot(lineTree, main="Line Mean clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
dev.off()

if(my.args$HCUT < Inf) {
  cat("Dropping lines clustering above", my.args$HCUT, "height\n")
  
  lineTreeSplit <- cutree(lineTree, h=my.args$HCUT)
  # Figure out which group to keep (the one with most lines)
  group.counts <- table(lineTreeSplit)
  keep.group <- names(group.counts)[which.max(group.counts)]
  keep.lines <- lineTreeSplit==keep.group
  
  cat("Dropping", sum(!keep.lines), "lines:", names(keep.lines)[!keep.lines], "\n")
  
  # IF THIS CHECK FAILS, JUST NEED TO REORDER keep.lines TO MATCH expr.table column IDs
  stopifnot(all(names(keep.lines)==colnames(expr.table)))
  
  expr.table <- expr.table[,keep.lines]
  
  cat(ncol(expr.table), "line means remain for further analysis.\n")
}


# --- Write output files --- #

cat("Writing remaining features for", my.args$MODE, "analysis to:", my.args$OUTDIR, "\n")

if(my.args$MODE == "mmc") {
  # The line mean output file should have same basename as EXPR input file, but with csv ending
  output.file <- basename(my.args$EXPR)
  output.file <- sub("[.]txt$", "", output.file)
  output.file <- paste0(my.args$OUTDIR, output.file, ".csv")
  cat("Writing filtered line means in CSV format to:", output.file, "\n")
  
  # Replace commas in feature names with "|" (impt for ambiguous features in transposon/microbiome analysis)
  row.names(expr.table) <- gsub(",", "|", row.names(expr.table))
  expr.table <- cbind(GENE=row.names(expr.table), expr.table)
  
  # Write out CSV file
  write.table(expr.table, output.file, sep=",", row.names=F, col.names=T, quote=F)
} else if(my.args$MODE == "eqtl") {
  # The pheno output file should have same basename as EXPR input file, but without _line_means suffix and .pheno instead of .txt
  pheno.file <- basename(my.args$EXPR)
  pheno.file <- sub("[.]txt$", "", pheno.file)
  pheno.file <- sub("[_]line[_][md][ei][af][nf]s$", "", pheno.file)
  pheno.file <- paste0(my.args$OUTDIR, pheno.file, ".pheno")
  
  cat("Outputting plink pheno format in", pheno.file, "\n")
  
  # Change column names to match "line_XXX" format for plink
  colnames(expr.table) <- sub("^X", "", colnames(expr.table))
  colnames(expr.table) <- sub("^DGRP[_]", "", colnames(expr.table))
  if(!grepl("^line[_]", colnames(expr.table)[1])) {
    colnames(expr.table) <- paste0("line_", colnames(expr.table))
  }
  
  # Strip out any trailing _F or _M parts of line columns
  # If this is Pooled line means but not SexAvg this will cause an error!
  colnames(expr.table) <- sub("[_][FM]$", "", colnames(expr.table))
  if(any(duplicated(colnames(expr.table)))) {
    stop("More than one column per line ID - if these are Pooled line means use corresponding SexAvg version instead!")
  }
  
  # Get rid of spaces and commas in row names (replace with underscores)
  row.names(expr.table) <- gsub(" ", "_", row.names(expr.table))
  row.names(expr.table) <- gsub(",", "_", row.names(expr.table))
  
  # Now transpose the matrix (that is the format plink wants)
  expr.table <- t(expr.table)
  
  # Append two columns for Family and Individual IDs (same thing in this case)
  expr.table <- cbind(FID=row.names(expr.table), IID=row.names(expr.table), as.data.frame(expr.table))
  write.table(expr.table, pheno.file, row.names=F, col.names=T, sep=" ", quote=F)
  cat("\n")
  
  if(my.args$PERM > 0) {
    for(perm in 1:my.args$PERM) {
      # Determine the output file path and name
      perm.path <- paste0(my.args$OUTDIR, "Perm", perm, "/")
      perm.file <- paste0(perm.path, basename(pheno.file))
      cat("Permuting into:", perm.file, "\n")
      
      # Make sure this subdir exists
      if(!file.exists(perm.path)) {
        dir.create(path=perm.path, recursive=T, showWarnings = F)
      }
      
      # Copy line means and shuffle the line IDs
      pheno.data <- expr.table
      line.shuf <- sample.int(n=nrow(pheno.data), replace=F)
      pheno.data$FID <- pheno.data$FID[line.shuf]
      pheno.data$IID <- pheno.data$IID[line.shuf]
      stopifnot(all(pheno.data$FID == pheno.data$IID))
      stopifnot(sum(duplicated(pheno.data$FID))==0)
      
      # Write out
      write.table(pheno.data, perm.file, col.names=T, row.names=F, quote=F, sep=" ")
    }
    cat("\n")
  }
  
} else {
  # WGCNA format - tab-delimited text, but transpose to sample x gene matrix
  # The output file should have same basename as EXPR input file, but with _wgcna suffix
  output.file <- basename(my.args$EXPR)
  output.file <- sub("[.]txt$", "", output.file)
  output.file <- paste0(my.args$OUTDIR, output.file, "_wgcna.txt")
  cat("Writing filtered line means in tab-delimited text format to:", output.file, "\n")
  
  # Transpose the expr.table to be Sample by Gene
  expr.table <- as.data.frame(t(expr.table))
  expr.table <- cbind(LINE=row.names(expr.table), expr.table)
  
  # Write out CSV file
  write.table(expr.table, output.file, sep="\t", row.names=F, col.names=T, quote=F)
}

cat("Script completed successfully!\n\n")
print(proc.time())
