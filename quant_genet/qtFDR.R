# ==========================
# = estimate empirical fdr =
# ==========================

args <- commandArgs(TRUE)

library(doMC)
registerDoMC(8)
fdr = 0.05

# phenotype file
# ============================================================

gene.list <- colnames(read.table("adjusted.line.means.pheno", header = TRUE, as.is = TRUE, row.names = 1)[, -1])

fdr.out <- foreach (i = 1:length(gene.list), .combine = rbind) %dopar% {
  
  # first read the observed p value
  obs.pval <- read.table(paste("gwas.", gene.list[i], ".qassoc", sep = ""),
                         header = TRUE, as.is = TRUE)
  res <- c(gene.list[i], fdr, NA, NA)
  
  if (nrow(obs.pval) > 0) {
    
    obs.pval <- obs.pval[order(obs.pval[, 9]), c(2, 9)]
    # read permuted data
    perm.pval <- data.frame(matrix(scan(con <- pipe(paste("cat perms/gwas.perm*.", gene.list[i], ".qassoc | awk \'$1 != \"CHR\" {print $2\"\\t\"$9}\'", sep = "")), what = ""), ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
    perm.pval[, 2] <- as.numeric(perm.pval[, 2])
    obs.pval <- cbind(obs.pval, rep(0, nrow(obs.pval)))
    for (j in 1:nrow(obs.pval)) {
      obs.pval[j, 3] <- sum(perm.pval[, 2] <= obs.pval[j, 2])/100/j
    }
    # identify the cutoff point
    thres <- which(cummax(obs.pval[, 3]) < fdr)
    if (length(thres) > 0) {
      res <- c(gene.list[i], obs.pval[max(thres), 3], paste(obs.pval[1:max(thres), 1], collapse = ","), paste(obs.pval[1:max(thres), 2], collapse = ","))
    } 
  }
  
  cat(gene.list[i], "\n")
  res
}

write.table(fdr.out, file = "qt.fdr.out", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

sessionInfo()
