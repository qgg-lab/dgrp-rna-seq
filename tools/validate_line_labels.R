#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 8/16/16
# 
# GOAL: Use the allele count data from each sample library to verify
#   the correct identity of each sample
#   
# USAGE: genotype_samples.R [OPTIONS] BATCH1 [BATCH2 ...]
#  --samples=FILE specifies the sample info file (Default: sample_info.txt)
#  --path=PATH specifies the top-level directory containing allele counts (Default: allele_counts/)
#  --tgeno=FILE specifies full path to tgeno file from DGRP (Default: /home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/dgrp2.tgeno)
#

options(stringsAsFactors=F)

# Take the list of batches to analyze from command line params
usageStr <- "USAGE:\nRscript genotype_samples.R [OPTIONS] BATCH1 [BATCH2 ...]\n  --samples=FILE specifies the sample master table file (Default: sample_info.txt)\n  --path=PATH specifies the top-level directory containing allele counts (Default: allele_counts/)"
myArgs <- commandArgs(trailingOnly=T)
# Interactive code testing: 
# myArgs <- c("--tgeno=~/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/dgrp2.tgeno", "batch_151216")
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
  cat("Setting --path=allele_counts/ by default.\n")
  myParams$path <- "allele_counts/"
}
if(!grepl("^/",myParams$path)) {
  myParams$path <- paste0(myParams$path, "/")
}

if(is.null(myParams$tgeno)) {
  cat("Setting --tgeno=/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/dgrp2.tgeno by default.\n")
  myParams$tgeno <- "/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/dgrp2.tgeno"
}

cat("\n")

# Load the tgeno file
if(!file.exists(myParams$tgeno)) {
  stop("Could not find required file ", myParams$tgeno)
}
tgeno.table <- read.table(myParams$tgeno, header=T, sep=" ", row.names=3)
# Drop all MNP and INDEL rows
tgeno.table <- tgeno.table[grepl("[_]SNP$", row.names(tgeno.table)),]
line.cols <- grep("^line[_]", colnames(tgeno.table), value=T)
line.numbers <- as.integer(sub("^line[_]","",line.cols))
names(line.numbers) <- line.cols
cat("Loaded genotype data for", nrow(tgeno.table), "SNPs across", length(line.cols), "lines from", myParams$tgeno, "\n\n")

# Load the sample master table
# TO DO: Would be nice to make this optional, 
# i.e. if the file or LINE col is missing, then just skip the actual verification check,
# but still report the best match for each sample
if(!file.exists(myParams$samples)) {
  stop("Could not find required file ", myParams$samples)
}
sample.table <- read.table(myParams$samples, header=T, sep="\t", row.names=1)
# Make sure it contains LINE column
if(!("LINE" %in% colnames(sample.table))) {
  stop(myParams$samples, " does not contain a LINE column, cannot verify correct line identities")
}
# Check that all LINE numbers are in genotype table
if(length(setdiff(sample.table$LINE, line.numbers))) {
  cat("!!!!!!!!!!\nWARNING:", myParams$samples, "contains some line numbers not in genotype table:", setdiff(sample.table$LINE, line.numbers), "\n!!!!!!!!!!\n\n")
}

for(batch in batches) {
  cat("\nAnalyzing all samples in batch:", batch, "\n\n")
  
  batch.file <- paste0(batch, ".txt")
  batch.path <- paste0(myParams$path, batch, "/")
  if(!file.exists(batch.file)) {
    stop("Could not find required batch file: ", batch.file)
  }
  batch.table <- read.table(batch.file, sep="\t")
  # Only use the first two columns. First column points to library file, second column is sample ID
  colnames(batch.table)[1:2] <- c("LIBRARY","SAMPLE")
  # Make sure all sample IDs in sample.table
  if(length(setdiff(batch.table$SAMPLE, row.names(sample.table))) > 0) {
    missing.samples <- setdiff(batch.table$SAMPLE, row.names(sample.table))
    stop(length(missing.samples), " sample names in ", batch.file, " are not found in ", myParams$samples, ": ", paste(missing.samples, collapse=" "))
  }
  
  # Start table to fill in top two hits and their scores
  # Also add the annotated line and error for that even if it's not the best or next
  batch.verify <- data.frame(FLAG=character(0), BEST.LINE=integer(0), BEST.ERR=numeric(0), NEXT.LINE=integer(0), NEXT.ERR=numeric(0), ANNOT.LINE=integer(0), ANNOT.ERR=numeric(0), RATIO.LINE=integer(0), RATIO.ERR=numeric(0))
  
  # Loop over each sample
  for(i in 1:nrow(batch.table)) {
    lib <- batch.table$LIB[i]
    sample <- batch.table$SAMPLE[i]
    line <- sample.table[sample,"LINE"]
    
    lib.file <- paste0(batch.path, lib, "_counts.txt")
    if(!file.exists(lib.file)) {
      cat("!!!!!!!!!!\nWARNING: Could not find",lib.file,"SKIPPING library for",lib,"(LINE ",line,")\n!!!!!!!!!!\n\n")
    } else {
      cat("Verifying genotype of", lib, "(Sample =", sample, ", Line =", line, ") from counts in", lib.file,"\n")
      
      # Load the allele_counts table
      allele.table <- read.table(lib.file, header=T, sep="\t")
      stopifnot(ncol(allele.table)==7)
      colnames(allele.table)[c(6,7)] <- c("FILTER","COUNTS")
      
      row.names(allele.table) <- paste(allele.table[,1], allele.table[,2], "SNP", sep="_")
      if(!all(row.names(allele.table) %in% row.names(tgeno.table))) {
        missing.snps <- setdiff(row.names(allele.table), row.names(tgeno.table))
        cat("\n!!!!!!!!!!\nWARNING: Dropping count data for", length(missing.SNPs),"SNPs missing from genotype table.\n!!!!!!!!!!\n\n")
        allele.table <- allele.table[row.names(allele.table) %in% row.names(tgeno.table),]
      }
      
      # Drop rows with 0 counts for this library
      allele.table <- allele.table[allele.table$COUNTS != "0,0",]
      
      # Subset tgeno table to what's in allele.table only
      tgeno.sub <- tgeno.table[row.names(allele.table),]
      # Make sure Ref allele matches across tables, warn if this is not the case
      if(!all(tgeno.sub$ref == allele.table$REF)) {
        cat("!!!!!!!!!!\nWARNING: REF alleles do not match for", sum(tgeno.sub$ref != allele.table$REF), "SNPs\n!!!!!!!!!!\n\n")
      }
      
      # Split COUNTS into MAJOR and MINOR
      allele.table[,"MAJ.COUNT"] <- as.integer(unlist(lapply(strsplit(allele.table$COUNTS, split=",", fixed=T), "[", 1)))
      allele.table[,"MIN.COUNT"] <- as.integer(unlist(lapply(strsplit(allele.table$COUNTS, split=",", fixed=T), "[", 2)))
      allele.table[,"DEPTH"] <- allele.table$MAJ.COUNT + allele.table$MIN.COUNT
      allele.table[,"MIN.FREQ"] <- allele.table$MIN.COUNT / allele.table$DEPTH
      
      # Drop SNPs that have perfect 50/50 frequency in the RNA-seq library
      # Warn if >5% rows are biallelic
      if(mean(allele.table$MAJ.COUNT == allele.table$MIN.COUNT) > 0.05) {
        cat("!!!!!!!!!!\nWARNING:", sum(allele.table$MAJ.COUNT == allele.table$MIN.COUNT),"of",nrow(allele.table), "of genotyped SNPs have strong heterozygosity.\n!!!!!!!!!!\n\n")
      }
      allele.table <- allele.table[allele.table$MAJ.COUNT != allele.table$MIN.COUNT,]
      tgeno.sub <- tgeno.sub[row.names(allele.table),]
      
      # TESTING: ONLY Use SNPs with at least 3 reads confirming genotype, and <= 25% of reads with minor genotype
      allele.table <- allele.table[(allele.table[,"MAJ.COUNT"] > 3) & (allele.table[,"MIN.FREQ"] <= 0.25),]
      tgeno.sub <- tgeno.sub[row.names(allele.table),]
      
      # Fill in CALL column that will be matched against line columns in tgeno table
      # For SNPs where MAJOR == REF, CALL = 0
      allele.table[allele.table$MAJOR==tgeno.sub$ref,"CALL"] <- 0
      # For SNPs where MAJOR == ALT, CALL = 2
      allele.table[allele.table$MAJOR==tgeno.sub$alt,"CALL"] <- 2
      # Everything else gets CALL = -1 (always erroneous)
      # This happens when MAJOR != REF or ALT
      allele.table[is.na(allele.table$CALL),"CALL"] <- -1
      # Warn if >5% of SNPs don't match the REF or ALT allele in known genotype
      if(mean(allele.table$CALL == -1) > 0.05) {
        cat("!!!!!!!!!!\nWARNING:", sum(allele.table$CALL == -1),"of",nrow(allele.table), "of genotyped SNPs have genotype not maching REF or ALT in DGRP table.\n")
        cat("Distribution of allele counts at these sites:\n")
        cat(min(allele.table[allele.table$CALL==-1,"DEPTH"]), boxplot.stats(allele.table[allele.table$CALL==-1,"DEPTH"])$stats, max(allele.table[allele.table$CALL==-1,"DEPTH"]), "\n")
        cat("Distribution of flags at these sites:\n")
        print(table(allele.table[allele.table$CALL==-1,"FILTER"]))
        cat("!!!!!!!!!!\n\n")
      }
      # Drop these SNPs since they won't contribute to match rate for any line
      allele.table <- allele.table[allele.table$CALL != -1,]
      tgeno.sub <- tgeno.sub[row.names(allele.table),]
      
      # WARN if much lower than typical amount, report total SNPs used
      if(nrow(allele.table) < (4*(10**5))) {
        cat("!!!!!!!!!!\nWARNING: < 500K SNPs were reliably genotyped in", lib, "\n!!!!!!!!!!\n")
      }
      cat(nrow(allele.table), "SNPs remain for verifying genotype of", lib, "\n")
      
      # Loop over each line, compute the error rate for all calls
      # Error rate = mismatched calls / total calls
      line.mean.error <- unlist(lapply(line.cols, function(j){
        # Which alleles in this line have genotype data, and are not segregating within line?
        line.geno <- tgeno.sub[,j] %in% c("0","2")
        # Should be able to use at least 80% of SNPs for each line
        if(mean(line.geno) < 0.8) {
          cat("!!!!!!!!!!\nWARNING:", j, "only has usable genotype data for", round(mean(line.geno)*100), "% of SNPS =", sum(line.geno), "SNPs for", j, "\n!!!!!!!!!!\n")
        }
        mean(allele.table[line.geno,"CALL"] != as.integer(tgeno.sub[line.geno,j]))
      }))
      names(line.mean.error) <- line.cols
      
      line.mean.error <- sort(line.mean.error, decreasing=F)
      
      best.line <- line.numbers[names(line.mean.error)[1]]
      best.err <- line.mean.error[1]
      next.line <- line.numbers[names(line.mean.error)[2]]
      next.err <- line.mean.error[2]
      
      # Alterate approach: Assume the line ID is correct
      # Loop over all OTHER lines, and check the genotypes at just the SNPs that would differentiate these two lines
      # Take the ratio of incorrect genotypes (other / current), which should always be > 1
      cur.line.col <- names(line.numbers)[line.numbers==line]
      other.line.cols <- names(line.numbers)[line.numbers != line]
      line.err.ratio <- unlist(lapply(other.line.cols, function(j){
        # Identify which SNPs have genotype for cur.line.col and j
        pair.seg <- (tgeno.sub[,cur.line.col] %in% c("0","2")) & (tgeno.sub[,j] %in% c("0","2")) & (tgeno.sub[,cur.line.col] != tgeno.sub[,j])
        # NOTE: This is often going to be an extremely small number of SNPs!
        # The only case to WARN on is when it's equal to 0 or < 5
        if(sum(pair.seg) == 0) {
          cat("!!!!!!!!!!\nWARNING: There are NO reliably genotyped SNPs known to segregate between", cur.line.col, "and", j, "\n!!!!!!!!!!\n")
          return(0)
        } else if(sum(pair.seg) < 5) {
          cat("!!!!!!!!!!\nWARNING: Only",sum(pair.seg),"reliably genotyped SNPs known to segregate between", cur.line.col, "and", j, "\n!!!!!!!!!!\n")
        }
        # Count the number that are incorrect for cur
        cur.err <- sum(allele.table[pair.seg,"CALL"] != as.integer(tgeno.sub[pair.seg,cur.line.col]))
        # Count the number that are incorrect for j
        j.err <- sum(allele.table[pair.seg,"CALL"] != as.integer(tgeno.sub[pair.seg,j]))
        # Return ratio of j.err : cur.err
        j.err / cur.err
      }))
      names(line.err.ratio) <- other.line.cols
      
      line.err.ratio <- sort(line.err.ratio, decreasing=F)
      best.ratio.line <- line.numbers[names(line.err.ratio)[1]]
      best.ratio <- line.err.ratio[1]
      
      # Check 1: best.line should be same as intended line
      if(best.line != line) {
        # TO DO: Make this a warning condition instead of an error condition if ratio method is more reliable
        cat("\n!!!!!!!!!!\nERROR:", lib, "is annotated as line", line, "but best match is line", best.line, "\n!!!!!!!!!!\n")
        batch.verify[lib,"FLAG"] <- "ERR"
      } else if((next.err / best.err) < 8) {
        # Check 2: next.err shoud be >8x higher than best.err
        cat("\n!!!!!!!!!!\nWARNING:", lib, "matches correct line", line, "but next best line", next.line, "has similar error rate.\n!!!!!!!!!!\n")
        batch.verify[lib,"FLAG"] <- "WARN"
      } else {
        cat("[OK]", lib, "matches line", best.line, "\n")
        batch.verify[lib,"FLAG"] <- "OK"
      }
      # Check 2: best.ratio should be > 1
      if(best.ratio <= 1) {
        cat("\n!!!!!!!!!!\nERROR:", lib, "is annotated as line", line, "but", best.ratio.line, "has error ratio =",best.ratio,"for segregating SNPs.\n!!!!!!!!!!\n")
        batch.verify[lib,"FLAG"] <- "ERR"
        # TO DO: COULD USE THIS TO UN-FLAG WARNINGS ABOVE IFF THE (best.ratio.line, line) = (best.line, next.line)
      }
      cat("\n")
      
      # Fill in table
      batch.verify[lib,"BEST.LINE"] <- best.line
      batch.verify[lib,"BEST.ERR"] <- best.err
      batch.verify[lib,"NEXT.LINE"] <- next.line
      batch.verify[lib,"NEXT.ERR"] <- next.err
      batch.verify[lib,"ANNOT.LINE"] <- line
      batch.verify[lib,"ANNOT.ERR"] <- line.mean.error[names(line.numbers)[line.numbers==line]]
      batch.verify[lib,"RATIO.LINE"] <- best.ratio.line
      batch.verify[lib,"RATIO.ERR"] <- best.ratio
    }
  }

  cat("\nSummary of batch", batch, ":\n")
  cat(sum(batch.verify$FLAG == "OK"), "libraries were verified as correct genotype.\n")
  if(any(batch.verify$FLAG == "WARN")) {
    cat(sum(batch.verify$FLAG == "WARN"), "libraries had next-best score warnings:\n")
    print(batch.verify[batch.verify$FLAG=="WARN",])
  }
  if(any(batch.verify$FLAG == "ERR")) {
    cat(sum(batch.verify$FLAG == "ERR"), "libraries had wrong genotype!\n")
    print(batch.verify[batch.verify$FLAG=="ERR",])
  }
  
  out.file <- paste0(batch.path, "genotype_check.txt")
  cat("Writing complete verification results for batch",batch,"to:",out.file,"\n")
  write.table(cbind(SAMPLE=row.names(batch.verify), batch.verify), out.file, sep="\t", row.names=F, quote=F)
  
  cat("\n\n")
}












