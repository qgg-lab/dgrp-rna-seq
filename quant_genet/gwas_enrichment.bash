# make permutation
# ============================================================

~/software/R-3.2.2/bin/Rscript -e 'set.seed(1); perm.seq <- matrix(ncol = 200, nrow = 100); for (i in 1:100) { perm.seq[i, ] <- sample(1:200); }; write.table(perm.seq, file = "perm.seq", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ");'

~/software/R-3.2.2/bin/Rscript -e 'perm.seq <- read.table("perm.seq", header = FALSE, as.is = TRUE); pheno <- read.table("adjusted.line.means.pheno", header = TRUE, as.is = TRUE); for (i in 1:100) { write.table(cbind(pheno[, 1:2], pheno[unlist(perm.seq[i, ]), -(1:2)]), file = paste("permPheno/qt.perm", i, ".pheno", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " "); cat(i, "\n"); }' 2>&1 &

# gwas, retain p < 1e-5
# ============================================================

~/software/plink-v1.90b3w/plink --silent --bfile /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/freeze2.200line.common --pheno adjusted.line.means.pheno --all-pheno --assoc --pfilter 0.00001 --out gwas > gwas.log 2>&1 

# permutation GWAS
# ============================================================

for i in `seq 1 100`
do
  srun ~/software/plink-v1.90b3w/plink --silent --bfile /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/freeze2.200line.common --pheno permPheno/qt.perm"$i".pheno --all-pheno --assoc --pfilter 0.00001 --out perms/gwas.perm"$i" > perms/gwas.perm"$i".log 2>&1 &
done

# estimate emprical FDR for the trait GWAS
# ============================================================

~/software/R-3.2.2/bin/Rscript qtFDR.R
