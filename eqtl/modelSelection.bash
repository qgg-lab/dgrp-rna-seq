# =======================================
# = run model selection to filter eQTLs =
# =======================================

# 1. extract SNPs from mapped eQTLs and get eQTL gene pairs
# ============================================================

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/ExprF.0.05.fdr.results.txt | sort -k1,1 | join -t $'\t' - <(sort ../female.genvar.txt) | grep -v NA | cut -f 3 | sed 's/,/\n/g' > modelSelect/female.eqtl.list

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/ExprM.0.05.fdr.results.txt | sort -k1,1 | join -t $'\t' - <(sort ../male.genvar.txt) | grep -v NA | cut -f 3 | sed 's/,/\n/g' > modelSelect/male.eqtl.list

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/ExprF.0.05.fdr.results.txt | sort -k1,1 | join -t $'\t' - <(sort ../female.genvar.txt) | grep -v NA > modelSelect/female.eqtl.pair

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/ExprM.0.05.fdr.results.txt | sort -k1,1 | join -t $'\t' - <(sort ../male.genvar.txt) | grep -v NA > modelSelect/male.eqtl.pair

# 2. extract genotypes from plink file
# ============================================================

~/software/plink-1.07-x86_64/plink --noweb --silent --bfile /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/freeze2.200line.common --extract modelSelect/female.eqtl.list --recode12 --transpose --out modelSelect/female.eqtl &

~/software/plink-1.07-x86_64/plink --noweb --silent --bfile /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/freeze2.200line.common --extract modelSelect/male.eqtl.list --recode12 --transpose --out modelSelect/male.eqtl &

# 3. get TSS information for all genes
#    this is used as a second criterion to select eQTLs
#    when association alone cannot distinguish two eQTLs
# ============================================================

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/combined-gene-info.txt | awk -F "\t" '{print $1"\t"$7"\t"$8}' | perl -wne 'chomp $_; @line = split /\t/, $_; @info = split /:/, $line[1]; @bound = split /-/, $info[1]; if ($line[2] eq "+") { print $line[0], "\t", $info[0], "\t", $line[2], "\t", $bound[0], "\n"; } else { print $line[0], "\t", $info[0], "\t", $line[2], "\t", $bound[1], "\n"; }' > modelSelect/gene.tss

# 4. phenotypes
# /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/combined_samples_known_novel_fpkm_VR_eQTL_F.pheno
# /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/combined_samples_known_novel_fpkm_VR_eQTL_M.pheno
# ============================================================

~/software/R-3.2.2/bin/Rscript modelSelect.R modelSelect/female.eqtl.tped modelSelect/female.eqtl.tfam /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/combined_samples_known_novel_fpkm_VR_eQTL_F.pheno modelSelect/female.eqtl.pair modelSelect/gene.tss 8 0.00001 20 modelSelect/female.modelSelect.RData > modelSelect/female.modelSelect.Rout 2>&1 &

~/software/R-3.2.2/bin/Rscript modelSelect.R modelSelect/male.eqtl.tped modelSelect/male.eqtl.tfam /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/plink/combined_samples_known_novel_fpkm_VR_eQTL_M.pheno modelSelect/male.eqtl.pair modelSelect/gene.tss 8 0.00001 20 modelSelect/male.modelSelect.RData > modelSelect/male.modelSelect.Rout 2>&1 &

# 5. output data and extract eQTL
# ============================================================

Rscript -e 'load("modelSelect/female.modelSelect.RData"); write.table(mc.select.eqtl, row.names = F, col.names = F, quote = F, sep = "\t")' | sed 's/VAR_//g' | sort -k1,1 | join -t $'\t' <(tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/combined-gene-info.txt | cut -f 1,7 | sort -k1,1) - > female.eqtl

Rscript -e 'load("modelSelect/male.modelSelect.RData"); write.table(mc.select.eqtl, row.names = F, col.names = F, quote = F, sep = "\t")' | sed 's/VAR_//g' | sort -k1,1 | join -t $'\t' <(tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/combined-gene-info.txt | cut -f 1,7 | sort -k1,1) - > male.eqtl

# 6. extract information for table S4
# ============================================================

Rscript -e 'load("modelSelect/female.modelSelect.RData"); write.table(mc.select.eqtl, row.names = F, col.names = F, quote = F, sep = "\t")' | sed 's/VAR_//g' | sort -k1,1 | join -t $'\t' <(tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/combined-gene-info.txt | cut -f 1,7,8 | sort -k1,1) - | join -t $'\t' - <(cut -f 1,3,4 modelSelect/female.eqtl.pair | sort -k1,1) | perl -wne 'chomp $_; @line = split /\t/, $_; ($chr, $pos) = split /:/, $line[1]; @bound = split /-/, $pos; @selected = split /,/, $line[3]; @eqtl_selected = (); for (my $i = 0; $i <= $#selected; $i++) { ($thischr, $thispos, $thistype) = split /_/, $selected[$i]; if ($thischr ne $chr) { push(@eqtl_selected, $selected[$i]."(trans-inter-arm)"); }  elsif ($thischr eq $chr && ( abs($thispos-$bound[0]) <= 1000 || abs($thispos - $bound[1]) <= 1000)) { push(@eqtl_selected, $selected[$i]."(cis)"); } else { push(@eqtl_selected, $selected[$i]."(trans-intra-arm)"); }  } @cis = (); @cisp = (); @intra = (); @intrap = (); @inter = (); @interp = (); @snps = split /,/, $line[4]; @pval = split /,/, $line[5]; for (my $i = 0; $i <= $#snps; $i++) { ($thischr, $thispos, $thistype) = split /_/, $snps[$i]; if ($thischr ne $chr) { push(@inter, $snps[$i]); push(@interp, $pval[$i]); }  elsif ($thischr eq $chr && ( abs($thispos-$bound[0]) <= 1000 || abs($thispos - $bound[1]) <= 1000)) { push(@cis, $snps[$i]); push(@cisp, $pval[$i]); } else { push(@intra, $snps[$i]); push(@intrap, $pval[$i]); }  } print $line[0], "\t", $chr, "\t", $bound[0], "\t", $bound[1], "\t", $line[2], "\t", join(",", @eqtl_selected), "\t", join(",", @cis), "\t", join(",", @cisp), "\t", join(",", @intra), "\t", join(",", @intrap), "\t", join(",", @inter), "\t", join(",", @interp), "\n";' > female.eqtl.select.info.txt 


Rscript -e 'load("modelSelect/male.modelSelect.RData"); write.table(mc.select.eqtl, row.names = F, col.names = F, quote = F, sep = "\t")' | sed 's/VAR_//g' | sort -k1,1 | join -t $'\t' <(tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/genes/combined-gene-info.txt | cut -f 1,7,8 | sort -k1,1) - | join -t $'\t' - <(cut -f 1,3,4 modelSelect/male.eqtl.pair | sort -k1,1) | perl -wne 'chomp $_; @line = split /\t/, $_; ($chr, $pos) = split /:/, $line[1]; @bound = split /-/, $pos; @selected = split /,/, $line[3]; @eqtl_selected = (); for (my $i = 0; $i <= $#selected; $i++) { ($thischr, $thispos, $thistype) = split /_/, $selected[$i]; if ($thischr ne $chr) { push(@eqtl_selected, $selected[$i]."(trans-inter-arm)"); }  elsif ($thischr eq $chr && ( abs($thispos-$bound[0]) <= 1000 || abs($thispos - $bound[1]) <= 1000)) { push(@eqtl_selected, $selected[$i]."(cis)"); } else { push(@eqtl_selected, $selected[$i]."(trans-intra-arm)"); }  } @cis = (); @cisp = (); @intra = (); @intrap = (); @inter = (); @interp = (); @snps = split /,/, $line[4]; @pval = split /,/, $line[5]; for (my $i = 0; $i <= $#snps; $i++) { ($thischr, $thispos, $thistype) = split /_/, $snps[$i]; if ($thischr ne $chr) { push(@inter, $snps[$i]); push(@interp, $pval[$i]); }  elsif ($thischr eq $chr && ( abs($thispos-$bound[0]) <= 1000 || abs($thispos - $bound[1]) <= 1000)) { push(@cis, $snps[$i]); push(@cisp, $pval[$i]); } else { push(@intra, $snps[$i]); push(@intrap, $pval[$i]); }  } print $line[0], "\t", $chr, "\t", $bound[0], "\t", $bound[1], "\t", $line[2], "\t", join(",", @eqtl_selected), "\t", join(",", @cis), "\t", join(",", @cisp), "\t", join(",", @intra), "\t", join(",", @intrap), "\t", join(",", @inter), "\t", join(",", @interp), "\n";' > male.eqtl.select.info.txt 

# counts
awk -F "\t" '$7 != ""' female.eqtl.select.info.txt | wc -l
awk -F "\t" '$7 != ""' male.eqtl.select.info.txt | wc -l

awk -F "\t" '$9 != "" || $11 != ""' female.eqtl.select.info.txt | wc -l
awk -F "\t" '$9 != "" || $11 != ""' male.eqtl.select.info.txt | wc -l

awk -F "\t" '$11 != ""' female.eqtl.select.info.txt | wc -l
awk -F "\t" '$11 != ""' male.eqtl.select.info.txt | wc -l

# 7. also permform model selection for TE expression
# ============================================================
tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/ActivityF.0.05.fdr.results.txt | awk '$3 != "NA"' | cut -f 3 | sed 's/,/\n/g' > modelSelect/female.te.eqtl.list

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/ActivityM.0.05.fdr.results.txt | awk '$3 != "NA"' | cut -f 3 | sed 's/,/\n/g' > modelSelect/male.te.eqtl.list

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/ActivityF.0.05.fdr.results.txt | awk '$3 != "NA"' > modelSelect/female.te.eqtl.pair

tail -n+2 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/ActivityM.0.05.fdr.results.txt | awk '$3 != "NA"' > modelSelect/male.te.eqtl.pair

~/software/plink-1.07-x86_64/plink --noweb --silent --bfile /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/freeze2.200line.common --extract modelSelect/female.te.eqtl.list --recode12 --transpose --out modelSelect/female.te.eqtl &

~/software/plink-1.07-x86_64/plink --noweb --silent --bfile /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/freeze2.200line.common --extract modelSelect/male.te.eqtl.list --recode12 --transpose --out modelSelect/male.te.eqtl &

# no tss for TE, use chr5
# ============================================================

cat <(head -n 1 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/combined_transposon_filtered_rpm_DNAAdj_eQTL_F.pheno | cut -d " " -f 3- | sed 's/ /\n/g') <(head -n 1 /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/combined_transposon_filtered_rpm_DNAAdj_eQTL_M.pheno | cut -d " " -f 3- | sed 's/ /\n/g') | sort | uniq | awk '{print $1"\t5\t+\t1"}' > modelSelect/te.tss

# select

~/software/R-3.2.2/bin/Rscript modelSelect.R modelSelect/female.te.eqtl.tped modelSelect/female.te.eqtl.tfam /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/combined_transposon_filtered_rpm_DNAAdj_eQTL_F.pheno modelSelect/female.te.eqtl.pair modelSelect/te.tss 8 0.00001 20 modelSelect/female.te.modelSelect.RData > modelSelect/female.te.modelSelect.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript modelSelect.R modelSelect/male.te.eqtl.tped modelSelect/male.te.eqtl.tfam /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Post/transposons/plink/combined_transposon_filtered_rpm_DNAAdj_eQTL_M.pheno modelSelect/male.te.eqtl.pair modelSelect/te.tss 8 0.00001 20 modelSelect/male.te.modelSelect.RData > modelSelect/male.te.modelSelect.Rout 2>&1 &

# test for overlap between gene eQTLs and TE eQTLs
# ============================================================

# there are a total of 1932427 variants
# in females, there are 66091 eQTLs in females and 104012 in males
# 2353 TE eQTLs in females and 4696 in males
# sort modelSelect/female.eqtl.list | uniq | wc -l
# sort modelSelect/male.eqtl.list | uniq | wc -l
# sort modelSelect/female.te.eqtl.list | uniq | wc -l
# sort modelSelect/male.te.eqtl.list | uniq | wc -l

# save data
# ============================================================

Rscript -e 'load("modelSelect/female.te.modelSelect.RData"); write.table(mc.select.eqtl, row.names = F, col.names = F, quote = F, sep = "\t")' | sed 's/VAR_//g'
Rscript -e 'load("modelSelect/male.te.modelSelect.RData"); write.table(mc.select.eqtl, row.names = F, col.names = F, quote = F, sep = "\t")' | sed 's/VAR_//g'


