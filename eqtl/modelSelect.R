## Perform model selection for eQTL
# core functions

# check if the addition of a genotype makes the genotype matrix funny,
# meaning that by removing missing genotypes,
# one geno becomes invariant
# or totally correlated with each other.

check.poly <- function(geno.matrix) {
  geno.matrix <- t(na.omit(t(geno.matrix)))
  if (ncol(geno.matrix) == 0) { # HAPPEN when no individuals remian
    return (FALSE)
  }
  alf <- rowMeans(geno.matrix)
  if (sum(alf == 0 | alf == 2) > 0) {
    return (FALSE)
  } else {
    if (nrow(geno.matrix) >= 3) {
      geno.matrix <- geno.matrix[-1, ]
      geno.corr <- cor(t(geno.matrix))
      if (sum(abs(geno.corr[upper.tri(geno.corr)]) >= 1-1e-10) > 0) {
        return (FALSE)
      } else {
        return (TRUE)
      }
    } else {
      return (TRUE)
    }
  }
}
  
select.eqtl <- function(snp.geno, exp.pheno, tss.chr, tss.pos, snp.pos, p.in, max.snp) {
  final.fit <- lm(exp.pheno ~ 1)
  # add a placeholder to the end of the genotype matrix
  # this ensures that the subsetting of matrix can work until the last snp

  snp.geno <- rbind(snp.geno, rep(0, length(exp.pheno)))

  # and the reason to set up model.geno is to make sure that
  # the model.geno matrix can be added
  # in check.poly(), the first row is always deleted
  # 1.5 is to make sure that it does not interfere with the condition
  # alf  == 0 | alf == 2
  model.geno <- rep(1.5, length(exp.pheno))
  rownames(snp.geno)[nrow(snp.geno)] = "placeholder"
  
  # initialize the loop
  n.snp <- 0
  
  while (n.snp < max.snp & nrow(snp.geno) > 0) {

    # intialize this search
    pre.p <- 1
    pre.geno <- rep(0, length(exp.pheno))
    best.fit <- NA
    pre.dist <- 30000000 # maximum possible length, bigger than any chromosome

    # loop through snps that are still in the snp.geno matrix
    for (i in 1:(nrow(snp.geno) - 1)) {
      snp.name <- rownames(snp.geno)[i]
      # get genotype
      eval(parse(text = paste(snp.name, " = as.integer(snp.geno[", i, ", ])", sep = "")))
      
      if (eval(parse(text = paste("check.poly(rbind(model.geno, ", snp.name, "))", sep = "")))) {
        eval(parse(text = paste("fit.add <- update(final.fit, . ~ . + ", snp.name, ")", sep = "")))
        if (is.na(coefficients(fit.add)[snp.name])) {
          add.p <- 1
        } else {
          add.p <- summary(fit.add)$coefficients[nrow(summary(fit.add)$coefficients), 4]
        }
        # check if they are perfectly correlated
        if (eval(parse(text = paste("sum(", snp.name, " != pre.geno, na.rm = T)", sep = ""))) == 0 | eval(parse(text = paste("sum(", snp.name, " != 2-pre.geno, na.rm = T)", sep = ""))) == 0) {
          add.p <- min(add.p, pre.p)
          pre.p <- add.p
        }
      
        snp.info <- snp.pos[i, ]
        if (snp.info[1] == tss.chr) {
          snp.dist <- abs(tss.pos - as.integer(snp.info[2]))
        } else {
          snp.dist <- 30000000
        }
        if (add.p < pre.p) {
          add.snp <- snp.name
          pre.p <- add.p
          pre.dist <- snp.dist
          pre.geno <- eval(parse(text = snp.name))
          best.fit <- fit.add
        } else {
          if (add.p == pre.p) {
            if (snp.dist < pre.dist) {
              add.snp <- snp.name
              pre.p <- add.p
              pre.dist <- snp.dist
              pre.geno <- eval(parse(text = snp.name))
              best.fit <- fit.add
            }
          }
        }
      }
    }

    # when there is one SNP that can be added
    if (pre.p < p.in) {
      final.fit <- best.fit
      model.geno <- rbind(model.geno, snp.geno[add.snp, ])
      snp.geno <- snp.geno[setdiff(rownames(snp.geno), add.snp), ]
      if (is.vector(snp.geno)) {
        n.snp <- max.snp + 1
        snp.geno <- rbind(snp.geno, snp.geno) # this is to make sure that nrow(snp.geno) still works in the next iter
      } else {
        n.snp <- n.snp + 1
      }
    } else {
      n.snp <- max.snp + 1
    }
  }
  return(final.fit)
}
      



args = commandArgs(TRUE)
# takes seven arguments
# 1. genotype file
# 2. tfam fle
# 3. phenotype file
# 4. eqtl pair to be selected
# 5. gene boundary
# args <- c("modelSelect/female.eqtl.tped", "modelSelect/female.eqtl.tfam", "/home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/combined_samples_known_novel_fpkm_VR_eQTL_F.pheno", "modelSelect/female.eqtl.pair", "modelSelect/gene.tss", 5, 0.00001, 20, "modelselect/out.RData")
# FBgn0000024


tped.file <- args[1]
tfam.file <- args[2]
pheno.file <- args[3]
eqtl.file <- args[4]
gene.tss.file <- args[5]
n.cpu <- as.numeric(args[6])
p.in <- as.numeric(args[7])
max.snp <- as.numeric(args[8])
out.file <- args[9]

library(doMC)
registerDoMC(n.cpu)

sessionInfo()

# read genotype
tfam <- read.table(tfam.file, header = FALSE, as.is = T)
n.line <- nrow(tfam)
tped <- matrix(scan(tped.file, what = ""), ncol = 4 + 2*n.line, byrow = TRUE)

snp.name <- paste("VAR_", tped[, 2], sep = "")
snp.loc <- tped[, c(1, 4)]
chr.name <- c("2L", "2R", "3L", "3R", "X", "4")
snp.loc[, 1] <- chr.name[as.numeric(snp.loc[, 1])]
rownames(snp.loc) <- snp.name

geno.code <- matrix(4 - (as.numeric(tped[, seq(from = 5, length = n.line, by = 2)]) + as.numeric(tped[, seq(from = 6, length = n.line, by = 2)])), ncol = n.line)
geno.code[geno.code == 4] <- NA
rownames(geno.code) <- snp.name

# read phenotype
pheno <- read.table(pheno.file, header = TRUE, as.is = TRUE, row.names = 1)
pheno <- pheno[tfam[, 1], ]

# read eqtl pair informaiton
eqtl.pair <- na.omit(read.table(eqtl.file, header = FALSE, as.is = TRUE, na.strings = "NA"))

# read.tss
gene.tss <- read.table(gene.tss.file, header = FALSE, as.is = TRUE, row.names = 1)

# model selection
mc.select.eqtl <- foreach(i = 1:nrow(eqtl.pair), .combine = rbind) %dopar% {
  snps <- paste("VAR_", unlist(strsplit(eqtl.pair[i, 3], split = ",")), sep = "")
  if (length(snps) >= 2) {
    snp.geno <- geno.code[snps, ]
    exp.pheno <- pheno[, gsub("-", ".", eqtl.pair[i, 1])]
    tss.chr <- gene.tss[eqtl.pair[i, 1], 1]
    tss.pos <- gene.tss[eqtl.pair[i, 1], 3]
    snp.pos <- snp.loc[snps, ]
    selected.model<- select.eqtl(snp.geno, exp.pheno, tss.chr, tss.pos, snp.pos, p.in, max.snp)
    model.terms <- names(coefficients(selected.model))
    result <- c(eqtl.pair[i, 1], paste(model.terms[2:length(model.terms)], collapse = ","))
  } else{
    result <- c(eqtl.pair[i, 1], snps)
  }
  result
}

save(mc.select.eqtl, file = out.file)
