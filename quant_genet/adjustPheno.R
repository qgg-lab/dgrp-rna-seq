# =================================================================
# = adjust line means for covariates (inversions, wolbachia, PCs) =
# =================================================================

args <- commandArgs(TRUE)
# args <- c("line.id.txt", "lineMeans", "adjustData.RData", "reportData/adjust.pheno.RData", "adjusted.line.means.pheno")
load(args[3])

# function to do adjustment
# ============================================================
 
adjustPheno <- function(raw.pheno) {

# raw.pheno is a two column data.frame where the first is
# line id and the second is raw phenotype

  common.line <- intersect(rownames(wolba), raw.pheno[, 1])
  rownames(raw.pheno) <- raw.pheno[, 1]
  # get data
  pheno.data <- cbind(raw.pheno[common.line, 2], wolba[common.line, 1], inv[common.line, c("In_2L_t", "In_2R_NS", "In_3R_P", "In_3R_K", "In_3R_Mo")], pcs[common.line, ])
  colnames(pheno.data)[2] <- "wolba"
  fit.form <- "pheno.data[, 1] ~ 1"
  if (length(unique(pheno.data[, "wolba"])) > 1) {
    fit.form <- paste(fit.form, " + factor(wolba)")
  }
  if (length(unique(pheno.data[, "In_2L_t"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_2L_t)")
  }
  if (length(unique(pheno.data[, "In_2R_NS"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_2R_NS)")
  }
  if (length(unique(pheno.data[, "In_3R_P"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_P)")
  }
  if (length(unique(pheno.data[, "In_3R_K"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_K)")
  }
  if (length(unique(pheno.data[, "In_3R_Mo"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_Mo)")
  }
  fit.form <- paste(fit.form, " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", sep = "")

  lm.fit <- lm(formula(fit.form), data = pheno.data)
  
  covar.coef <- summary(lm.fit)$coefficients
  
  pheno.adjust <- residuals(lm.fit) + covar.coef[1]
  
  return(pheno.adjust)
  
}

# main codes to process phenotypes
# ============================================================

# profile the directory
# ============================================================

line.id <- scan(file = args[1], what = "", quiet = TRUE)

female.trait.files <- list.files(path = args[2], pattern = "*\\.female\\.csv")
male.trait.files <- list.files(path = args[2], pattern = "*\\.male\\.csv")

female.traits <- gsub("\\.female\\.csv", "", female.trait.files)
male.traits <- gsub("\\.male\\.csv", "", male.trait.files)

# loop through traits to load data
# ============================================================

female.line.means <- matrix(nrow = length(line.id), ncol = length(female.traits))
for (i in 1:length(female.traits)) {
  trait.name <- female.traits[i]
  # read data
  female.line.means[, i] <- read.csv(paste(args[2], "/", trait.name, ".female.csv", sep = ""), header = FALSE, as.is = TRUE, row.names = 1)[line.id, 1]
}
colnames(female.line.means) <- female.traits

male.line.means <- matrix(nrow = length(line.id), ncol = length(male.traits))
for (i in 1:length(male.traits)) {
  trait.name <- male.traits[i]
  # read data
  male.line.means[, i] <- read.csv(paste(args[2], "/", trait.name, ".male.csv", sep = ""), header = FALSE, as.is = TRUE, row.names = 1)[line.id, 1]  
}
colnames(male.line.means) <- male.traits

# make adjustment
# ============================================================

female.adjust.pheno <- female.line.means

for (i in 1:ncol(female.line.means)) {
  
  this.trait.result <- adjustPheno(data.frame(line <- as.character(line.id), pheno <- female.line.means[, i], stringsAsFactors = F))
  female.adjust.pheno[, i] <- this.trait.result[line.id]
  
}

colnames(female.adjust.pheno) <- paste("female.", female.traits, sep = "")

male.adjust.pheno <- male.line.means

for (i in 1:ncol(male.line.means)) {
  
  this.trait.result <- adjustPheno(data.frame(line <- as.character(line.id), pheno <- male.line.means[, i], stringsAsFactors = F))
  male.adjust.pheno[, i] <- this.trait.result[line.id]
  
}
colnames(male.adjust.pheno) <- paste("male.", male.traits, sep = "")

# write data
# ============================================================

write.table(rbind(c("FID", "IID", colnames(female.adjust.pheno), colnames(male.adjust.pheno)), cbind(line.id, line.id, female.adjust.pheno, male.adjust.pheno)), file = args[5], sep = " ", col.names = F, row.names = F, quote = F)
