#!/bin/bash
#
# LJE - 5/25/16
# All analysis steps for applying the pipeline to known_novel_genes/combined_samples_known_novel_counts.txt
# TO DO: Generalize this pipeline so it can be invoked on any directory/count table

# SESSION SETTINGS
USENODE="-x node1"

# ASSEMBLY splits the analysis by which build of novel transcriptome to use

# [DEPRECATED]
# rep = reproducible novel genes only)
# ASSEMBLY="rep"
# STUB=cuffcompare.reproducible.novel

# Main version of the analysis
# all = using all novel genes found by cufflinks
ASSEMBLY="all"
STUB=cuffcompare.${ASSEMBLY}.novel
GTFFILE=$STUB.gtf
INFOFILE=${STUB}-gene-info.txt
COUNTFILE=combined_samples_${ASSEMBLY}_novel_counts.txt
MYPATH="genes"
COUNTSTUB="combined_samples_known_novel"


# ONLY NEED TO RUN THIS ONCE:
# Variant density of the KNOWN genes
sbatch $USENODE -o compute_variant_density_known.Rout compute_variant_density.R ~/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-transcriptome-r5.57-plus-r6.11-backport.gff known_gene

mkdir -p $MYPATH

# Compute variant density tables for novel genes
sbatch -o compute_variant_density_novel_${ASSEMBLY}.Rout compute_variant_density.R $GTFFILE

# First combine the known and novel feature counts (also variant density tables):
Rscript combine_known_novel_gene_counts.R $COUNTFILE $STUB $MYPATH &> $MYPATH/combine_known_novel_gene_counts.Rout

# Copy the combined GTF and related info files to Resources directory
RESPATH="/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/"
mkdir -p $RESPATH
cp known_all_novel_genes/combined_gene_models.gtf $RESPATH
cp known_all_novel_genes/combined-gene-info.txt $RESPATH
cp known_all_novel_genes/combined_gene_variant_rates.txt $RESPATH
cp known_all_novel_genes/combined_samples_known_novel_counts.txt $RESPATH


# TO DO: Try VOOM here as well

# Standard normalization (now includes TMM!)
# Normalizes to column sums
# TO DO: Keep LOW and RARE features in at this stage? 
# Do normalized output files still get FLAG column?
LOGFILE=$MYPATH/normalize_expression.Rout
normalize_expression.R COUNTS=$MYPATH/$COUNTSTUB"_counts.txt" OUTDIR=$MYPATH GENES=$MYPATH/combined-gene-info.txt CUTOFF=FIT &> $LOGFILE

# Copy the normalized expression data to Resources directory
cp known_all_novel_genes/combined_samples_known_novel_fpkm.txt $RESPATH


# --- Genetic Variance Model --- #

mkdir -p $MYPATH/genVar

# Primary H^2 analysis (Within Sex) - the main analysis uses FPKM values, skips NQ transform, corrects for variant rate AND Wolbachia
# NOTE: Parallelization auto-detects the number of CPUs given to jobs
sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_VR_WolAdj.Rout gen_var_model.R EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_variant_rates.txt TAG=VR OUTDIR=$MYPATH/genVar/

# DEPRECATED
# # Run permutations of primary analysis to make sure H2 is not technical artifact
# # TO DO: NEED TO RE-RUN THIS!
# # SKIPPING FOR UPDATED SCRIPT FOR NOW
# mkdir -p $MYPATH/genVar/Perm/
# for i in {1..10}
# do
#  sbatch $USENODE -c 20 -o $MYPATH/genVar/Perm/gen_var_model_fpkm_VR_WolAdj_Perm${i}.Rout gen_var_model.R EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_variant_rates.txt TAG=VR OUTDIR=$MYPATH/genVar/Perm/ PERM=$i
# done

# Pooled Sex H^2 Analysis
# RE-RUNNING WITH FLAG VERSION ON 1/18/17
# TO DO: Should do VR correction within sex? Or include a SEX:LINEREG ixn term?
sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_VR_Pooled_Wol.Rout gen_var_model.R POOLED=TRUE EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_variant_rates.txt TAG=VR OUTDIR=$MYPATH/genVar/

# eQTL H^2 Analysis (Within Sex) - Same as primary analysis but includes additional model terms for major inversions and cryptic relatedness PCs
# RE-RUNNING WITH FLAG VERSION ON 1/18/17
# These line means should be used for eQTL mapping
sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_VR_eQTL.Rout gen_var_model.R EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW EQTL=T LINEREG=$MYPATH/combined_gene_variant_rates.txt TAG=VR OUTDIR=$MYPATH/genVar/

# Alternate strategies for comparison - these can be run in parallel with the command above

# DEPRECATED
# # Try subsets of VR counts (poly and indel only)...
# sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_polyVR_WolAdj.Rout gen_var_model.R EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_poly_variant_rates.txt TAG=polyVR OUTDIR=$MYPATH/genVar/
# sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_indelVR_WolAdj.Rout gen_var_model.R EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_indel_variant_rates.txt TAG=indelVR OUTDIR=$MYPATH/genVar/
# # This must be run AFTER running the VR-corrected genetic variance model
# # TO DO: UPDATE THIS SCRIPT
# sbatch $USENODE -o gene_variants_vs_expr.Rout gene_variants_vs_expr.R
# 
# # TO DO: If this correction is mild enough in both cases, it should become part of the standard build
# # There should be one No-VR correction version just for comparison, but otherwise all downstream processing should include the VR correction

# TO DO: Rename the log file from this, give it an extra stub like noVR?
# Try skipping the Variant Rate correction
sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_noVR_WolAdj.Rout gen_var_model.R EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW TAG=noVR OUTDIR=$MYPATH/genVar/

# Try skipping the Wolbachia correction
sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_VR_noWol.Rout gen_var_model.R WOLADJ=FALSE EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_variant_rates.txt TAG=VR OUTDIR=$MYPATH/genVar/

# DEPRECATED
# # Try RLE normalization of FPKM values
# sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_rle_VR_WolAdj.Rout gen_var_model.R EXPR=$MYPATH/$COUNTSTUB"_rle.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_variant_rates.txt TAG=VR OUTDIR=$MYPATH/genVar/
#
# # Try NQ normalization of FPKM values
# sbatch $USENODE -c 20 -o $MYPATH/genVar/gen_var_model_fpkm_NQ_VR_WolAdj.Rout gen_var_model.R NQ=TRUE EXPR=$MYPATH/$COUNTSTUB"_fpkm.txt" FILTER=LOW LINEREG=$MYPATH/combined_gene_variant_rates.txt TAG=VR OUTDIR=$MYPATH/genVar/

# Run comparison of main H2 analysis vs all alternate versions:
# TO DO: The first 4 comparisons should all be done on transposon and microbiome features too
# TO DO: Need a way to compare results from pooled vs single analysis
Rscript compare_H2_tables.R REF=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" ALT=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_NoWol_model_results.txt" &> $MYPATH/genVar/compare_H2_Wol_effects.txt
Rscript compare_H2_tables.R REF=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" ALT=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_eQTL_model_results.txt" &> $MYPATH/genVar/compare_H2_InvPC_effects.txt
# NOTE: ~500 genes are no longer significant in eQTL model, but ~100 genes become significant in this model - eQTL analysis should be run on genes significant in BOTH only!
# Rscript compare_H2_tables.R REF=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" ALT=$MYPATH/genVar/$COUNTSTUB"_rle_VR_WolAdj_model_results.txt" &> $MYPATH/genVar/compare_H2_rle_effects.txt
# Rscript compare_H2_tables.R REF=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" ALT=$MYPATH/genVar/$COUNTSTUB"_fpkm_NQ_VR_WolAdj_model_results.txt" &> $MYPATH/genVar/compare_H2_NQ_effects.txt
# Compare VR results vs no VR results
Rscript compare_H2_tables.R REF=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" ALT=$MYPATH/genVar/$COUNTSTUB"_fpkm_noVR_WolAdj_model_results.txt" &> $MYPATH/genVar/compare_H2_VR_effects.txt
# Compare full VR results vs the other subsets
# Rscript compare_H2_tables.R REF=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" ALT=$MYPATH/genVar/$COUNTSTUB"_fpkm_indelVR_WolAdj_model_results.txt" &> $MYPATH/genVar/compare_H2_indelVR_effects.txt
# Rscript compare_H2_tables.R REF=$MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" ALT=$MYPATH/genVar/$COUNTSTUB"_fpkm_polyVR_WolAdj_model_results.txt" &> $MYPATH/genVar/compare_H2_polyVR_effects.txt

# Run comparison of these model results to what I got before fixing issue with line mean extraction
# All model_results.txt tables should be the same
# Line Means for WolAdj only models should be very similar (most correlations > 0.9)
# Line means for eQTL models should be less correlated to what we had before
# NOTE: This also runs the comparison for microbiome and transposon features
sbatch $USENODE -o compare_new_genVar.Rout compare_new_genVar.R

# TO DO: Update the plotting script to split out known and novel features separately?
# For now, just use grep to split the appropriate files and run the plotting script separately?
# NOTE: The plotting script just ignores the FLAG column

# Plots from primary model
grep '^[FG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" > $MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_model_results.txt
grep '^[XG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_model_results.txt" > $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_model_results.txt

Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_model_results.txt &> $MYPATH/genVar/plot_all_heritability_FPKM.Rout
Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_model_results.txt &> $MYPATH/genVar/plot_known_heritability_FPKM.Rout
Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_model_results.txt &> $MYPATH/genVar/plot_novel_heritability_FPKM.Rout

# Plots from unadjusted model (for explanatory figures)
grep '^[FG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_noVR_WolAdj_model_results.txt" > $MYPATH/genVar/combined_samples_split_known_fpkm_noVR_WolAdj_model_results.txt
Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_split_known_fpkm_noVR_WolAdj_model_results.txt &> $MYPATH/genVar/plot_known_noVR_heritability_FPKM.Rout
# Run comparison on just the known genes
Rscript compare_H2_tables.R REF=$MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_model_results.txt ALT=$MYPATH/genVar/combined_samples_split_known_fpkm_noVR_WolAdj_model_results.txt &> $MYPATH/genVar/compare_H2_VR_effects_known_only.txt


# Plots from pooled model
grep '^[FG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_Pooled_model_results.txt" > $MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_Pooled_model_results.txt
grep '^[XG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_Pooled_model_results.txt" > $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_Pooled_model_results.txt

Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_known_novel_fpkm_VR_WolAdj_Pooled_model_results.txt &> $MYPATH/genVar/plot_all_pooled_FPKM.Rout
Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_Pooled_model_results.txt &> $MYPATH/genVar/plot_known_pooled_FPKM.Rout
Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_Pooled_model_results.txt &> $MYPATH/genVar/plot_novel_pooled_FPKM.Rout

# Plots from eQTL model
grep '^[FG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_eQTL_model_results.txt" > $MYPATH/genVar/combined_samples_split_known_fpkm_VR_eQTL_model_results.txt
grep '^[XG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_eQTL_model_results.txt" > $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_eQTL_model_results.txt

Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_known_novel_fpkm_VR_eQTL_model_results.txt &> $MYPATH/genVar/plot_all_eQTL_FPKM.Rout
Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_split_known_fpkm_VR_eQTL_model_results.txt &> $MYPATH/genVar/plot_known_eQTL_FPKM.Rout
Rscript plot_bicolor_hist.R $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_eQTL_model_results.txt &> $MYPATH/genVar/plot_novel_eQTL_FPKM.Rout

# Copy primary and eQTL model results to Resources directory
# (UPDATED COPIES MADE ON 2/25/18)
cp $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_F_line_means.txt" $RESPATH
cp $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_M_line_means.txt" $RESPATH
cp $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_eQTL_F_line_means.txt" $RESPATH
cp $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_eQTL_M_line_means.txt" $RESPATH

# U R HERE - 2/26/18
# Reran the remaining models above with updated genvar script
# Log reports look good, comparison script looks mostly good 
# (a few genes have low corr even in WolAdj only model, but most are >0.9)
# NEXT: proceed with downstream analysis


# --- Principal Component Analysis --- #

mkdir -p $MYPATH/pca
sbatch $USENODE -c 20 -o $MYPATH/pca/expression_pca_fpkm_VR_WolAdj_F.Rout expression_PCA.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt
sbatch $USENODE -c 20 -o $MYPATH/pca/expression_pca_fpkm_VR_WolAdj_M.Rout expression_PCA.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt


# --- Trait Correlation --- #

# TESTING WITH NEW FLAG COLUMN ON 1/18/17 - DOES NOT PASS THRU, CAN IT IGNORE?
# TO DO: This script drops the FLAG column as a "missing line", should handle this more explicitly
# Compute trait correlations using FPKM line means (both sexes, +/- Log2 transform for now)
# NOTE: Line means are now kept on log2 scale, so there is no need to use LOG2 option
# TO DO: Add option to do correlation across both sexes (correlate avg and diff?)
mkdir -p $MYPATH/traitCorr
sbatch $USENODE -o $MYPATH/traitCorr/trait_correlations_fpkm_Wol_F.Rout trait_correlations.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt
sbatch $USENODE -o $MYPATH/traitCorr/trait_correlations_fpkm_Wol_M.Rout trait_correlations.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt

# Also split the line means for known and novel genes, and try trait correlation analysis separately on each
grep '^[FG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_F_line_means.txt" > $MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_F_line_means.txt
grep '^[FG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_M_line_means.txt" > $MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_M_line_means.txt
grep '^[XG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_F_line_means.txt" > $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_F_line_means.txt
grep '^[XG]' $MYPATH/genVar/$COUNTSTUB"_fpkm_VR_WolAdj_M_line_means.txt" > $MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_M_line_means.txt

sbatch $USENODE -o $MYPATH/traitCorr/trait_correlations_split_known_fpkm_Wol_F.Rout trait_correlations.R EXPR=$MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_F_line_means.txt
sbatch $USENODE -o $MYPATH/traitCorr/trait_correlations_split_known_fpkm_Wol_M.Rout trait_correlations.R EXPR=$MYPATH/genVar/combined_samples_split_known_fpkm_VR_WolAdj_M_line_means.txt
sbatch $USENODE -o $MYPATH/traitCorr/trait_correlations_split_novel_fpkm_Wol_F.Rout trait_correlations.R EXPR=$MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_F_line_means.txt
sbatch $USENODE -o $MYPATH/traitCorr/trait_correlations_split_novel_fpkm_Wol_M.Rout trait_correlations.R EXPR=$MYPATH/genVar/combined_samples_split_novel_fpkm_VR_WolAdj_M_line_means.txt

# Plot all trait correlation distributions
Rscript plot_bicolor_hist.R $MYPATH/traitCorr/${COUNTSTUB}_fpkm_VR_WolAdj_F_traitCorr.txt &> $MYPATH/traitCorr/plot_correlations_FPKM_F.Rout
Rscript plot_bicolor_hist.R $MYPATH/traitCorr/${COUNTSTUB}_fpkm_VR_WolAdj_M_traitCorr.txt &> $MYPATH/traitCorr/plot_correlations_FPKM_M.Rout

Rscript plot_bicolor_hist.R $MYPATH/traitCorr/combined_samples_split_known_fpkm_VR_WolAdj_F_traitCorr.txt &> $MYPATH/traitCorr/plot_correlations_known_FPKM_F.Rout
Rscript plot_bicolor_hist.R $MYPATH/traitCorr/combined_samples_split_known_fpkm_VR_WolAdj_M_traitCorr.txt &> $MYPATH/traitCorr/plot_correlations_known_FPKM_M.Rout
Rscript plot_bicolor_hist.R $MYPATH/traitCorr/combined_samples_split_novel_fpkm_VR_WolAdj_F_traitCorr.txt &> $MYPATH/traitCorr/plot_correlations_novel_FPKM_F.Rout
Rscript plot_bicolor_hist.R $MYPATH/traitCorr/combined_samples_split_novel_fpkm_VR_WolAdj_M_traitCorr.txt &> $MYPATH/traitCorr/plot_correlations_novel_FPKM_M.Rout

# -- New Version -- #

# This uses a more standardized script and does linear regression
# TO DO: Add spearman correlations into this version?
# TO DO: Add plots to the core script here?
# TO DO: Lots of other potential improvements, see TO DO items in the script itself

mkdir -p $MYPATH/traitreg
# Run against individual gene features
sbatch $USENODE -o $MYPATH/traitreg/trait_regression_${COUNTSTUB}_fpkm_VR_WolAdj_F.Rout trait_regression.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt
sbatch $USENODE -o $MYPATH/traitreg/trait_regression_${COUNTSTUB}_fpkm_VR_WolAdj_M.Rout trait_regression.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt

# Run against PCs
sbatch $USENODE -o $MYPATH/traitreg/trait_regression_${COUNTSTUB}_fpkm_VR_WolAdj_F_PCs.Rout trait_regression.R EXPR=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_F_PCs.txt
sbatch $USENODE -o $MYPATH/traitreg/trait_regression_${COUNTSTUB}_fpkm_VR_WolAdj_M_PCs.Rout trait_regression.R EXPR=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_M_PCs.txt


# --- Run MMC --- #

# --- DEPRECATED - USE WGCNA INSTEAD --- #
#
# mkdir $MYPATH/mmc
# 
# # Seem to have had some success with Female gene expr but only after VERY stringent filtering
# # Need to explore this further and consider other methods
# 
# # Change this parameter depending on desired level of resolution
# # When testing MMC:
# RESOLUTION="--fast"
# # For final results, re-run with:
# RESOLUTION=""
# 
# # VAR=0.05 seems to work MUCH BETTER for F
# # TO DO: Try slightly higher VAR cutoff for M? (VAR=0.1 worked but maybe too strict)
# # ALSO NOTE: Before filtering line means, PC1.Perc.Var is highly correlated with specificity for gonads (Ovary in F, Testis in M) and NEGATIVELY correlated with variance
# # There is also a clear correlation between the larger tightly correlated clusters and low variance
# #  - my hope is that once we remove lower variance genes, we can get more equitable clustering
# 
# # Filter out LOW/RARE genes, low VAR genes, and filter avg expression a bit more stringently
# # NOTE: Threshold on M may need to be a bit higher to get similar # of genes
# filter_line_means.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt FILTER=RARE,LOW AVGEXPR=0 VAR=0.05 &> $MYPATH/mmc/filter_all_F.Rout
# filter_line_means.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt FILTER=RARE,LOW AVGEXPR=0 VAR=0.05 &> $MYPATH/mmc/filter_all_M.Rout
# 
# # Run MMC - Using Pearson Correlation (default)
# sbatch $USENODE --exclusive -o $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_F_mmc.log run_mmc.sh $RESOLUTION --correlation=pearson $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.csv
# sbatch $USENODE --exclusive -o $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_M_mmc.log run_mmc.sh $RESOLUTION --correlation=pearson $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.csv
# 
# # Run MMC - Using Spearman Correlation
# sbatch $USENODE --exclusive -o $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_F_mmc_scc.log run_mmc.sh $RESOLUTION --correlation=spearman --suffix=scc $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.csv
# sbatch $USENODE --exclusive -o $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_M_mmc_scc.log run_mmc.sh $RESOLUTION --correlation=spearman --suffix=scc $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.csv
# 
# # Plot results
# plot_mmc_results.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt MMC=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_mmc.csv &> $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_F_plot_mmc_results.Rout &
# plot_mmc_results.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt MMC=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_mmc.csv &> $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_M_plot_mmc_results.Rout &
# 
# plot_mmc_results.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt MMC=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_mmc_scc.csv &> $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_F_plot_mmc_scc_results.Rout &
# plot_mmc_results.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt MMC=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_mmc_scc.csv &> $MYPATH/mmc/${COUNTSTUB}_fpkm_Wol_M_plot_mmc_scc_results.Rout &
# 
# # Look at stats for clustering results vs various other gene properties
# sbatch -o $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_mmc_covariates.Rout analyze_covariates.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt CLUST=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_mmc.csv PCA=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_F_PC_gene_results.txt
# sbatch -o $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_mmc_covariates.Rout analyze_covariates.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt CLUST=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_mmc.csv PCA=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_M_PC_gene_results.txt
# 
# sbatch -o $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_mmc_scc_covariates.Rout analyze_covariates.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt CLUST=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_F_mmc_scc.csv PCA=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_F_PC_gene_results.txt
# sbatch -o $MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_mmc_scc_covariates.Rout analyze_covariates.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt CLUST=$MYPATH/mmc/${COUNTSTUB}_fpkm_VR_WolAdj_M_mmc_scc.csv PCA=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_M_PC_gene_results.txt
# 
# ---------------------------------- #


# --- Run WCGNA --- #

# This is an alternate Clustering Method to MMC
mkdir -p $MYPATH/wgcna

# Filter out LOW/RARE genes, and filter by avg expr, H2, and VAR as well:
filter_line_means.R MODE=wgcna EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt FILTER=RARE,LOW AVGEXPR=0 VAR=0.05 &> $MYPATH/wgcna/filter_all_F.Rout
filter_line_means.R MODE=wgcna EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt FILTER=RARE,LOW AVGEXPR=0 VAR=0.05 &> $MYPATH/wgcna/filter_all_M.Rout

# Run WGCNA on F and M
# TO DO: Plenty of room for improvement or alternate things to try here,
# See TO DO list on run_wgcna.R
sbatch $USENODE -c 20 -o $MYPATH/wgcna/run_wgcna_F.Rout run_wgcna.R EXPR=$MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means_wgcna.txt
sbatch $USENODE -c 20 -o $MYPATH/wgcna/run_wgcna_M.Rout run_wgcna.R EXPR=$MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means_wgcna.txt

# Plot results
nohup plot_mmc_results.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt MMC=$MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna.csv &> $MYPATH/wgcna/${COUNTSTUB}_fpkm_Wol_F_plot_wgcna_results.Rout &
nohup plot_mmc_results.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt MMC=$MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna.csv &> $MYPATH/wgcna/${COUNTSTUB}_fpkm_Wol_M_plot_wgcna_results.Rout &

sbatch -o $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna_covariates.Rout analyze_covariates.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_F_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt CLUST=$MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna.csv PCA=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_F_PC_gene_results.txt
sbatch -o $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna_covariates.Rout analyze_covariates.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_M_line_means.txt H2FILE=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_WolAdj_model_results.txt CLUST=$MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna.csv PCA=$MYPATH/pca/${COUNTSTUB}_fpkm_VR_WolAdj_M_PC_gene_results.txt

# Extract the list of genes for each cluster to a separate directory
MAXK=`awk -F',' '{print $2}' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna.csv | uniq | sort -nr | head -n 1`
for ((k=1;k<=MAXK;k++))
do
  # Extract the list of KNOWN gene IDs into a new wgcna subdir for each cluster
  mkdir -p $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_cluster${k}
  awk -F',' '($2 == '$k') {print $1}' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna.csv > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_cluster${k}/${COUNTSTUB}_fpkm_VR_WolAdj_F_cluster${k}_genes_all.txt
  grep '^FBgn' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_cluster${k}/${COUNTSTUB}_fpkm_VR_WolAdj_F_cluster${k}_genes_all.txt > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_cluster${k}/${COUNTSTUB}_fpkm_VR_WolAdj_F_cluster${k}_genes_known.txt
done

MAXK=`awk -F',' '{print $2}' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna.csv | uniq | sort -nr | head -n 1`
for ((k=1;k<=MAXK;k++))
do
  # Extract the list of KNOWN gene IDs into a new wgcna subdir for each cluster
  mkdir -p $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_cluster${k}
  awk -F',' '($2 == '$k') {print $1}' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna.csv > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_cluster${k}/${COUNTSTUB}_fpkm_VR_WolAdj_M_cluster${k}_genes_all.txt
  grep '^FBgn' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_cluster${k}/${COUNTSTUB}_fpkm_VR_WolAdj_M_cluster${k}_genes_all.txt > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_cluster${k}/${COUNTSTUB}_fpkm_VR_WolAdj_M_cluster${k}_genes_known.txt
done

# Background set for enrichment analysis:
awk -F',' '($1 != "Gene") {print $1}' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna.csv > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna_input_genes_all.txt
grep '^FBgn' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna_input_genes_all.txt > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_F_wgcna_input_genes_known.txt

awk -F',' '($1 != "Gene") {print $1}' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna.csv > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna_input_genes_all.txt
grep '^FBgn' $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna_input_genes_all.txt > $MYPATH/wgcna/${COUNTSTUB}_fpkm_VR_WolAdj_M_wgcna_input_genes_known.txt

# U R HERE - 2/28/18
# Ran all steps in this section and committed those changes
# TO DO: Compare figures to what we had before
# TO DO: Upload cluster and input gene lists to FlyMine and download GO, domains, pathway, publications, bdgp results
# I may have written a description of what I did previously - check scrapbook
# Also check uncommited deleted files in genes/wgcna/ for naming conventions that are compatible with my script for table S6


# --- eQTL Analysis (Plink) --- #

# This only needs to be run ONCE for all analyses (microbiome, transposons, genes, etc.)
# ALREADY RUN ON MICROBIOME
# setup_plink.sh
# sbatch $USENODE -c 20 -o build_SNP_gene_map.Rout build_SNP_gene_map.R SNP=freeze2.200line.common.snp 

# Running eQTL mapping separately on each sex, using the DNAAdj line means (activity)
# Could also run on expression with NO DNAAdj (might capture QTL influencing integration activity as well, but would be harder to interpret)
# Pooled model shows widespread sex by line effects, so no reason to use AvgSex here

# TO DO: Could split this into 2 scripts
# The first runs ./setup_plink.sh steps IF the required files are missing, filters line means, queues up the jobs
# The second runs the summarization and cleanup steps (could also confirm that everything ran correctly?)

# Reformat line means for plink input, and create permuted versions
# (starting with 10 permutations, should expand to 100 or 1000 later?)
PNUM=100
mkdir -p $MYPATH/plink
filter_line_means.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_eQTL_F_line_means.txt FILTER=RARE,LOW MODE=eqtl PERM=$PNUM &> $MYPATH/plink/plink_setup_F.Rout
filter_line_means.R EXPR=$MYPATH/genVar/${COUNTSTUB}_fpkm_VR_eQTL_M_line_means.txt FILTER=RARE,LOW MODE=eqtl PERM=$PNUM &> $MYPATH/plink/plink_setup_M.Rout

# Run the main eQTL mapping and permutations (can be run in parallel)
sbatch $USENODE -o $MYPATH/plink/plink_run_F.log run_plink.sh --genotypes=freeze2.200line.common --phenotypes=$MYPATH/plink/${COUNTSTUB}_fpkm_VR_eQTL_F.pheno --output=$MYPATH/plink/ExprF
sbatch $USENODE -o $MYPATH/plink/plink_run_M.log run_plink.sh --genotypes=freeze2.200line.common --phenotypes=$MYPATH/plink/${COUNTSTUB}_fpkm_VR_eQTL_M.pheno --output=$MYPATH/plink/ExprM
for ((perm=1;perm<=PNUM;perm++))
do
  sbatch $USENODE -o $MYPATH/plink/Perm${perm}/plink_run_F.log run_plink.sh --genotypes=freeze2.200line.common --phenotypes=$MYPATH/plink/Perm${perm}/${COUNTSTUB}_fpkm_VR_eQTL_F.pheno --output=$MYPATH/plink/Perm${perm}/ExprF
  sbatch $USENODE -o $MYPATH/plink/Perm${perm}/plink_run_M.log run_plink.sh --genotypes=freeze2.200line.common --phenotypes=$MYPATH/plink/Perm${perm}/${COUNTSTUB}_fpkm_VR_eQTL_M.pheno --output=$MYPATH/plink/Perm${perm}/ExprM
done

# THIS WAS EXPERIMENTAL, BUT DID NOT WORK - IGNORE FOR NOW
# # Run the eQTL mapping with LD pruning activated
# # (If this gets used, it will be applied AFTER FDR filter, so it doesn't to be run on all permutations)
# # Putting the output in a separate subdir just in case, but in general the output file names
# # should not conflict with those from primary run
# # UPDATE: THIS DID NOT WORK
# # Plink design is basically idiotic for this, as far as I can tell
# # It seems that the prune step needs to be run *independently* on EVERY SINGLE qassoc file
# # So we basically would need to re-run the entire eQTL stage, then run prune step on EVERY DAMN QASSOC FILE
# # No. Thanks.
# mkdir $MYPATH/plink/pruned
# sbatch $USENODE -o $MYPATH/plink/pruned/plink_run_F_pruned.log run_plink.sh --genotypes=freeze2.200line.common --phenotypes=$MYPATH/plink/${COUNTSTUB}_fpkm_VR_eQTL_F.pheno --output=$MYPATH/plink/pruned/ExprF --prune
# sbatch $USENODE -o $MYPATH/plink/pruned/plink_run_M_pruned.log run_plink.sh --genotypes=freeze2.200line.common --phenotypes=$MYPATH/plink/${COUNTSTUB}_fpkm_VR_eQTL_M.pheno --output=$MYPATH/plink/pruned/ExprM --prune
# # TO DO: Should be able to run the same tabulation script used below on the pruned version
# # Can then try another filtering step that only keeps the 5% FDR SNPs that also made it past pruning step

# TEMP CODE TO TABULATE INDIVIDUAL PERM DIRECTORIES
#
# First check that all are done
# tail Perm5?/*.log | more
# 
# for pdir in Perm5?
# do
#  echo $pdir
#  cd $pdir
#  sbatch -o eQTL_tabulate.Rout eQTL_tabulate.R ./
#  cd ..
# done
#
# # Check this task finished OK
# more Perm5?/eQTL_tabulate.Rout
# wc -l Perm5?/*.results.txt
# ls -lh Perm5?/*.results.txt
#
# # Remove the qassoc files
# for pdir in Perm6? Perm7? Perm8? Perm9? Perm100
# do
#  echo $pdir
#  cd $pdir
#  rm *.qassoc
#  cd ..
# done

# AFTER ALL PLINK JOBS FINISH:

# Combine all the qassoc files in each dir into a single file
# TO DO: Add options to eQTL_tabulate.R and eQTL_perm_fdr to only process a specific model group (can get better parallelization that way)
sbatch $USENODE -o $MYPATH/plink/eQTL_tabulate.Rout eQTL_tabulate.R $MYPATH/plink/

# Run the summarization script that builds table of all eQTLs passing 5% FDR
# TO DO: Add option or additional script to cleanup all the intermediate files
# One thing worth saving is the collection of ALL permutation indexes (table with a column showing the permuted line order for each permutation)
sbatch $USENODE -o $MYPATH/plink/eQTL_perm_fdr.Rout eQTL_perm_fdr.R $MYPATH/plink/

# Repeat the eQTL permutation testing at more stringent thresholds
sbatch $USENODE -o $MYPATH/plink/eQTL_perm_fdr.0.01.Rout eQTL_perm_fdr.R $MYPATH/plink/ 0.01
sbatch $USENODE -o $MYPATH/plink/eQTL_perm_fdr.0.005.Rout eQTL_perm_fdr.R $MYPATH/plink/ 0.005
sbatch $USENODE -o $MYPATH/plink/eQTL_perm_fdr.0.001.Rout eQTL_perm_fdr.R $MYPATH/plink/ 0.001

# Everything from here proceeds with 0.05 FDR, but may want to consider adapting to other FDR

# Clean-up the qassoc files at this point:
cd $MYPATH/plink
rm *.qassoc
for pdir in Perm*
do
  echo $pdir
  cd $pdir
  rm *.qassoc
  cd ..
done


# Draw genome-wide map of eQTLs (all, not just network links)
# TO DO: This script should have a required parameter to specify the path to the FDR result table
# TO DO: Should also incorporate this script into DGRPseq toolkit
# (Move this further down if adding in coloration by eQTL type or network link)
sbatch $USENODE -o $MYPATH/plink/draw_gene_eQTL_map.Rout draw_gene_eQTL_map.R

# Parse these eQTL results, tabulate cis vs trans eQTLs
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_summary_F.Rout eQTL_summary.R EQTL=$MYPATH/plink/ExprF.0.05.fdr.results.txt MAP=freeze2.200line.common.snp.gene.map
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_summary_M.Rout eQTL_summary.R EQTL=$MYPATH/plink/ExprM.0.05.fdr.results.txt MAP=freeze2.200line.common.snp.gene.map
# Repeat using 1% FDR
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_summary_F.0.01.fdr.Rout eQTL_summary.R EQTL=$MYPATH/plink/ExprF.0.01.fdr.results.txt MAP=freeze2.200line.common.snp.gene.map
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_summary_M.0.01.fdr.Rout eQTL_summary.R EQTL=$MYPATH/plink/ExprM.0.01.fdr.results.txt MAP=freeze2.200line.common.snp.gene.map

# Compute networks
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_network_F.Rout eQTL_network.R EQTL=$MYPATH/plink/ExprF.0.05.fdr.trans.eqtls.txt
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_network_M.Rout eQTL_network.R EQTL=$MYPATH/plink/ExprM.0.05.fdr.trans.eqtls.txt
# Repeat these with 1% FDR as well
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_network_F.0.01.fdr.Rout eQTL_network.R EQTL=$MYPATH/plink/ExprF.0.01.fdr.trans.eqtls.txt
sbatch -c 20 $USENODE -o $MYPATH/plink/eQTL_network_M.0.01.fdr.Rout eQTL_network.R EQTL=$MYPATH/plink/ExprM.0.01.fdr.trans.eqtls.txt

# Compare new and old networks
# (Changing this to show how much changed after correcting the issue with gen_var_model.R)
compare_eQTL_networks.R $MYPATH/plink/ExprF.0.05.fdr.cis.trans.network.txt OLD_GEN_VAR/known_all_novel_genes/plink/ExprF.0.05.fdr.cis.trans.network.txt &> $MYPATH/plink/compare.new.vs.old.ExprF.0.05.cis.trans.network.Rout
compare_eQTL_networks.R $MYPATH/plink/ExprM.0.05.fdr.cis.trans.network.txt OLD_GEN_VAR/known_all_novel_genes/plink/ExprM.0.05.fdr.cis.trans.network.txt &> $MYPATH/plink/compare.new.vs.old.ExprM.0.05.cis.trans.network.Rout

# Compare Female and Male networks (5% FDR)
compare_eQTL_networks.R $MYPATH/plink/ExprF.0.05.fdr.cis.trans.network.txt $MYPATH/plink/ExprM.0.05.fdr.cis.trans.network.txt &> $MYPATH/plink/compare.ExprF.vs.ExprM.0.05.cis.trans.network.Rout

# Compare 5% vs 1% FDR networks
compare_eQTL_networks.R $MYPATH/plink/ExprF.0.05.fdr.cis.trans.network.txt $MYPATH/plink/ExprF.0.01.fdr.cis.trans.network.txt &> $MYPATH/plink/compare.ExprF.0.05.vs.0.01.cis.trans.network.Rout
compare_eQTL_networks.R $MYPATH/plink/ExprM.0.05.fdr.cis.trans.network.txt $MYPATH/plink/ExprM.0.01.fdr.cis.trans.network.txt &> $MYPATH/plink/compare.ExprM.0.05.vs.0.01.cis.trans.network.Rout

# Compare Female and Male networks (1% FDR)
compare_eQTL_networks.R $MYPATH/plink/ExprF.0.01.fdr.cis.trans.network.txt $MYPATH/plink/ExprM.0.01.fdr.cis.trans.network.txt &> $MYPATH/plink/compare.ExprF.vs.ExprM.0.01.cis.trans.network.Rout


# --- Summary Scripts --- #

# Create most of the tables for S5
sbatch $USENODE -o $MYPATH/table_S5.Rout table_S5.R
# NOTE: Only GEO tables 3 and 4 changed here - 5D was just a precision change on one number

sbatch $USENODE -o $MYPATH/wgcna/table_S6_WGCNA.Rout table_S6_WGCNA.R

# Summarization of gene counts in different categories
# The first one did not change with genVar updates
./known_vs_novel_expr.R &> known_vs_novel_expr.Rout
# This second one changed the png output files, but they are probably not meaningful changes:
./array_vs_RNAseq.R &> array_vs_RNAseq.Rout
# DEPRECATED: This script compared the Rep vs All versions of novel transcriptome, but the Rep version has been removed
# ./novel_rep_vs_all.R &> novel_rep_vs_all.Rout

# Final summarization script
# Depends on eQTL results!
# Re-ran after genVar fix - only the eQTL numbers change
./gene_result_summary.R &> gene_result_summary.Rout

# U R HERE - 3/9/18
# Re-ran enrichment analysis for WGCNA clusters and rebuilt Supplemental tables
# TO DO: Update the Excel files for Trudy
# TO DO: Check for other scripts that were missing here
# 		One thing that's missing is the known vs novel correlation analysis

# Interesting note: after adding more permutations to eQTL analysis
# The number of expression features with at least 1 eQTL reduced, and overall number of eQTLs reduced
# But the size of the eQTL cis->trans network actually increased a bit
# Suggests eQTLs in the network are more reliable

# TO DO LIST FOR GENE ANALYSIS:
# 1) Check for extent of sex effects and sex by line interactions - we don't actually know if within-sex analysis is justified for these features
# 2) PCA to check for overall pattern
#
