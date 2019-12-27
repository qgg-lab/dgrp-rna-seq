#!/bin/bash
#
# LJE - 2/9/16
#
# Complete list of top-level commands for DGRP Baseline RNA-seq Alignment
# This directory should contain only the scripts/files necessary for the main alignment pipeline
# Small batch testing commands should go in DGRP_Baseline_RNAseq_Sandbox
#


# --- SESSION SETUP --- #

# Project-level parameters
# RUN THIS COMMAND WHENEVER STARTING A NEW SESSION
source project_info.sh

# NODE USAGE - CHANGE THIS TO ADAPT TO HYPERION LOAD
# In general, it's always good practice to exclude at least one node
# to keep it open for other users
# For example, exclude nodes 1 and 2 (use 3 and 4 only)
USENODE="-x node1"

# Define the batches to run alignment tasks on:
# By default, run on all batches defined in project_info.sh
# But if you are running the analysis piecemeal, e.g. if you already aligned some of your
# data but then sequenced additional samples to process, then you can define a subset of
# your samples through this variable.  All steps below that can be re-run for just the
# new batches will use RUNBATCHES instead of ALLBATCHES
RUNBATCHES=${ALLBATCHES[@]}

# Batch for downstream processing of merged libraries from the same sample
MERGED_BATCH="merged_samples"


# Info about each batch
# To run a step on a specific batch, run one of these BATCH=... commands
# Then run the inner part of the loop for that step

# BATCH LIST:
#
# - Batch #1 (Run on 10/31/14)
# One flowcell, all 8 lanes
# BATCH="batch_141031"
#
# - Batch #2 (Run on 3/18/15)
# One flowcell, only 5 lanes
# BATCH="batch_150318"
#
# - Batch #3 (Run on 4/27/15)
# Full flowcell, all 8 lanes
# BATCH="batch_150427"
#
# - Batch #4 (Run on 5/27/15)
# Full flowcell, all 8 lanes
# BATCH="batch_150527"
#
# - Batch #5 (Run on 6/13/15)
# !!! MASKING THIS BATCH !!!
# ALL RELATED FILES MOVED TO BAD_BATCHES/
# There were major quality issues with this run
# BATCH="batch_150613"
#
# - Batch #6 (Run on 7/7/15)
# This is a double batch - two flow-cells that were run together on same date
# For technical artifact checking, this should be considered as batch6A and batch6B based on the flow cells
# BATCH="batch_150707"
#
# - Batch #7 (Run on 9/17/15)
# Another double batch, should be split into batch6A and batch6B for technical artifact checking
# BATCH="batch_150917"
#


# PRE-CHECK: Make sure all raw files in a batch are present
# NOTE: Raw Fastq files should be gzipped already!
for BATCH in ${RUNBATCHES[@]}
do
  ./run_batch_task.sh --batch=$BATCH --task='ls -l $RAW_FASTQ_PATH/$sampleID.fastq.gz' | more
done


# --- STEP 1: Run fastqc on raw input sequences (OPTIONAL) --- #

for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p fastqc/$BATCH
  logfile="fastqc/$BATCH/run_fastqc_slurm_%A_%a.out"
  ./run_batch_task.sh --array --threads=1 --batch=$BATCH --sbatchopt="-o $logfile $USENODE" --task='fastqc --extract --outdir fastqc/$argBatch $RAW_FASTQ_PATH/$sampleID.fastq.gz'
done

# Summarize results
./fastqc_summary.sh fastqc


# --- STEP 2: Run cutadapt (REQUIRED) --- #

# This step trims out adapter sequences that occur when insert < read length

for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $TRIMMED_FASTQ_PATH/$BATCH
  logfile="$TRIMMED_FASTQ_PATH/$BATCH/run_cutadapt_slurm_%A_%a.log"
  ./run_batch_task.sh --array --threads=1 --batch=$BATCH --sbatchopt="-o $logfile $USENODE" --task='./run_cutadapt.sh --output=$TRIMMED_FASTQ_PATH/$argBatch --sample=$sampleID'
done

# Summarize results
for BATCH in ${RUNBATCHES[@]}
do
  ./summarize_cutadapt.sh $TRIMMED_FASTQ_PATH/$BATCH/run_cutadapt_slurm_*.log > $TRIMMED_FASTQ_PATH/$BATCH.cutadapt.summary.txt
done

# Collect QC stats for trimmed fastq.gz files
for BATCH in ${RUNBATCHES[@]}
do
  outfile="$TRIMMED_FASTQ_PATH/$BATCH/fastq.gz.stats"
  errfile="$TRIMMED_FASTQ_PATH/$BATCH/fastq.gz.stats.err"
  sbatch -o $outfile -e $errfile ./file_stats.sh $TRIMMED_FASTQ_PATH/$BATCH/*.fastq.gz
done

cat $TRIMMED_FASTQ_PATH/batch_*/fastq.gz.stats > $TRIMMED_FASTQ_PATH/fastq.gz.stats


# --- STEP 3: Filter Ribosomal reads --- #

# Align trimmed reads to rRNA sequences using BWA
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $BWA_RIBO_PATH/$BATCH
  logfile="$BWA_RIBO_PATH/$BATCH/run_bwa_batch_ribo_slurm_%A_%a.out"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE -c $DEFAULT_THREAD_COUNT" --task='./run_bwa.sh --sample=$sampleID --batch=$argBatch'
done

# Collect QC stats for ribofiltered fastq.gz files:
for BATCH in ${RUNBATCHES[@]}
do
  outfile="$RIBO_FILTERED_PATH/$BATCH/fastq.gz.stats"
  errfile="$RIBO_FILTERED_PATH/$BATCH/fastq.gz.stats.err"
  sbatch -o $outfile -e $errfile ./file_stats.sh $RIBO_FILTERED_PATH/$BATCH/*.fastq.gz
done

cat $RIBO_FILTERED_PATH/batch_*/fastq.gz.stats > $RIBO_FILTERED_PATH/fastq.gz.stats

# Summarize how many reads aligned to each rRNA sequence
for BATCH in ${RUNBATCHES[@]}
do
  logfile="$BWA_RIBO_PATH/$BATCH/run_bwa_count_slurm_%A_%a.out"
  ./run_batch_task.sh --array --threads=1 --batch=$BATCH --sbatchopt="-o $logfile $USENODE" --task='bwa_count.py --verbose --nosplit bwa_rRNA/$argBatch/$sampleID/output.bam > bwa_rRNA/$argBatch/$sampleID/target_counts.txt'
done


# --- STEP 4: Filter Microbial Reads --- #

# Filter microbial reads using fast BWA alignment against DGRP_Microbiome database
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $BWA_MICROBE_PATH/$BATCH
  logfile="$BWA_MICROBE_PATH/$BATCH/run_bwa_microbe_slurm_%A_%a.out"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE -c $DEFAULT_THREAD_COUNT" --task='./run_bwa.sh --sample=$sampleID --batch=$argBatch --mode=microbe'
done

# Collect QC stats for microbefiltered fastq.gz files:
for BATCH in ${RUNBATCHES[@]}
do
  outfile="$MICROBE_FILTERED_PATH/$BATCH/fastq.gz.stats"
  errfile="$MICROBE_FILTERED_PATH/$BATCH/fastq.gz.stats.err"
  sbatch -o $outfile -e $errfile ./file_stats.sh $MICROBE_FILTERED_PATH/$BATCH/*.fastq.gz
done

cat $MICROBE_FILTERED_PATH/batch_*/fastq.gz.stats > $MICROBE_FILTERED_PATH/fastq.gz.stats

# Count reads aligning to each microbe sequence
for BATCH in ${RUNBATCHES[@]}
do
  logfile="$BWA_MICROBE_PATH/$BATCH/run_bwa_count_slurm_%A_%a.out"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE" --task='bwa_count.py --nosplit --verbose bwa_microbe/$argBatch/$sampleID/output.bam > bwa_microbe/$argBatch/$sampleID/target_counts.txt'
done

# --- ADDED 10/17/16 --- #
# Re-align microbe-aligning reads to Drosophila genome with STAR to filter out ambiguous reads
for BATCH in ${RUNBATCHES[@]}
do
  logfile="$BWA_MICROBE_PATH/$BATCH/run_star_realign_microbe_slurm_%A_%a.out"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE -c $DEFAULT_THREAD_COUNT" \
   --task='run_star.sh --sample=$sampleID --input=bwa_microbe/$argBatch/$sampleID/ --output=bwa_microbe/$argBatch/$sampleID/ --suffix="_mapped"'
done

# Could probably just delete these but keeping them for now just in case
gzip $BWA_MICROBE_PATH/$BATCH/*/Unmapped.out.mate1

# Extract IDs that align to the real chromosomes, not U(Extra), as that may contain microbial contamination
for BATCH in ${RUNBATCHES[@]}
do
  logfile="$BWA_MICROBE_PATH/$BATCH/run_bwa_count_star_realign_slurm_%A_%a.out"
  run_batch_task.sh --array --threads=1 --batch=$BATCH --sbatchopt="-o $logfile $USENODE" \
   --task='./extract_microbe_star_IDs.sh $BWA_MICROBE_PATH/$argBatch/$sampleID'
done

# Count reads aligning to each microbe sequence but NOT Drosophila genome
for BATCH in ${RUNBATCHES[@]}
do
  logfile="$BWA_MICROBE_PATH/$BATCH/run_bwa_count_star_realign_masked_slurm_%A_%a.out"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE" --task='bwa_count.py --nosplit --verbose --mask=bwa_microbe/$argBatch/$sampleID/reads_aligned_Dmel_star.txt bwa_microbe/$argBatch/$sampleID/output.bam > bwa_microbe/$argBatch/$sampleID/target_filtered_counts.txt'
done

# U R HERE - 10/18/16

# TO DO: See how much filtered and original counts differ at organism, genus level, etc. (do as part of microbe_counts.R overhaul)
#		- Also to change in microbe_counts.R - combine sample-level counts into libraries BEFORE condensing to organism, species, etc.
#		- Test replicate library correlation at assembly level
#		- Have options to control which version of counts are being used
#		- Maybe a separate script to compare filtered and original counts at all levels (also identify if some samples have more filtering than others?)
#


# --- STEP 5: Filter Repeat Element Reads --- #

# BWA Alignment of reads against RepBase:
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $BWA_REPEAT_PATH/$BATCH
  logfile="$BWA_REPEAT_PATH/$BATCH/run_bwa_repeat_slurm_%A_%a.out"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE -c $DEFAULT_THREAD_COUNT" --task='./run_bwa.sh --sample=$sampleID --batch=$argBatch --mode=repeat'
done

for BATCH in ${RUNBATCHES[@]}
do
  outfile="$REPEAT_FILTERED_PATH/$BATCH/fastq.gz.stats"
  errfile="$REPEAT_FILTERED_PATH/$BATCH/fastq.gz.stats.err"
  sbatch -o $outfile -e $errfile $USENODE ./file_stats.sh $REPEAT_FILTERED_PATH/$BATCH/*.fastq.gz
done

cat $REPEAT_FILTERED_PATH/batch_*/fastq.gz.stats > $REPEAT_FILTERED_PATH/fastq.gz.stats

# Count reads aligning to each repeat sequence
for BATCH in ${RUNBATCHES[@]}
do
  logfile="$BWA_REPEAT_PATH/$BATCH/run_bwa_count_slurm_%A_%a.out"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE" --task='bwa_count.py --verbose --nosplit bwa_repeat/$argBatch/$sampleID/output.bam > bwa_repeat/$argBatch/$sampleID/target_counts.txt'
done

# Run fastqc on the final reads meant for genomic alignment
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p fastqc_repfiltered/$BATCH
  logfile="fastqc_repfiltered/$BATCH/run_fastqc_slurm_%A_%a.out"
  ./run_batch_task.sh\
   --array\
   --threads=1\
   --batch=$BATCH\
   --sbatchopt="-o $logfile $USENODE"\
   --task='fastqc --extract --outdir fastqc_repfiltered/$argBatch $REPEAT_FILTERED_PATH/$argBatch/${sampleID}_filtered.fastq.gz'
done

./fastqc_summary.sh fastqc_repfiltered


# --- STEP 6: Align to Reference Genome using STAR --- #

# Align the fully filtered reads to reference genome using STAR
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $STAR_PATH/$BATCH
  logfile="$STAR_PATH/$BATCH/run_star_batch_slurm_%A_%a.log"
  ./run_batch_task.sh\
   --array\
   --batch=$BATCH\
   --sbatchopt="-o $logfile $USENODE -c $DEFAULT_THREAD_COUNT"\
   --task='run_star.sh --sample=$sampleID --batch=$argBatch --lane=$lane --name=$sampleName --flowcell=$flowcell'
done

# Compress the Unmapped.out.mate1 files
for BATCH in ${RUNBATCHES[@]}
do
  nohup gzip $STAR_PATH/$BATCH/*/Unmapped.out.mate1 &
done


# --- STEP 7: Count Reads in Known Genes using HTSeq-Count --- #

# HTSeq-count: Assemble individual sample read counts in all known gene features
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $HTSEQ_COUNT_PATH/$BATCH
  logfile="$HTSEQ_COUNT_PATH/$BATCH/run_htseq_batch_star_slurm_%A_%a.log"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="$USENODE -o $logfile -c 1" --task='run_htseq.sh --sample=$sampleID --batch=$argBatch --suffix=STAR_counts'
done


# --- STEP 8: QC/Validation --- #

# Run various QC steps at this stage
# WARNING: There are still some steps in these scripts that map batch names to simpler Batch1, 2, etc. numbering scheme
# 			These are specific to this project and need to be generalized or removed!

# NOTE: Not using DGRPseq versions of these scripts b/c there are too many project-specific checks and corrections being made here
# See Leips and DGRP_3WK projects for standardized versions of these...

# Collect sample stats
sbatch -o collect_sample_QC_stats.Rout ./collect_sample_QC_stats.R ${ALLBATCHES[@]}

# Generate plots of sample stats
sbatch -o plot_sample_QC_stats.Rout ./plot_sample_QC_stats.R

# Count observed alleles at all known SNPs
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $PILEUP_PATH/$BATCH/
  logfile="$PILEUP_PATH/$BATCH/run_count_alleles_slurm_%A_%a.log"
  run_batch_task.sh --array --batch=$BATCH --sbatchopt="-o $logfile $USENODE" --task='count_alleles_rna.sh --sample=$sampleID --outpath='$PILEUP_PATH/$BATCH/' --header '$STAR_PATH/$BATCH/'$sampleID/Aligned.sortedByCoord.out.bam'
done

# Run the genotyping script individually on each batch
for BATCH in ${RUNBATCHES[@]}
do
  logfile="$PILEUP_PATH/$BATCH/genotype_samples.Rout"
  sbatch -o $logfile $USENODE genotype_samples.R $BATCH
done

# Draw genotype error plots
./genotype_validation_plots.R ${ALLBATCHES[@]} &> genotype_validation_plots.Rout

# Sex validation script using PCA and LDA
validate_sex_labels.R ${ALLBATCHES[@]} &> validate_sex_labels.Rout

# Build master sample table for downstream processing
# THIS SCRIPT MUST RUN TO COMPLETION BEFORE STARTING NEXT SCRIPTS
sbatch -o build_sample_table.Rout ./build_sample_table.R ${ALLBATCHES[@]}
# NOTE: ALL samples still >=5M depth, majority (88%) still >= 10M


# --- STEP 9: COMBINE LIBRARIES FOR REPLICATED SAMPLES --- #

# Aggregate gene counts from HTSeq-count, repbase sequences, and microbiomes by library and sample
# These can all be run in parallel
# For main gene feature count table, now using the shared toolset version!
sbatch -o build_count_table.Rout ~/Tools/DGRPseq/build_count_table.R
sbatch -o build_lib_count_table.Rout ~/Tools/DGRPseq/build_lib_count_table.R
sbatch -o repeat_counts.Rout repeat_counts.R
sbatch -o microbe_counts.Rout -c 20 microbe_counts.R
sbatch -o microbe_filtered_counts.Rout -c 20 microbe_counts.R COUNTS=target_filtered_counts.txt OUTSTUB=filtered_counts.txt
Rscript microbe_count_compare.R &> microbe_count_compare.Rout

# Merge unmapped fastq files for Trinity analysis and aligned BAM files for cufflinks
sbatch $USENODE -o merge_unaligned_fastq.log merge_unaligned_fastq.sh
sbatch $USENODE -o merge_star_aligned_bam.log merge_star_aligned_bam.sh

# Make merged_samples batch file
grep -v '^SAMPLE[[:space:]]LIB' sample_master_table.txt | awk -F'\t' '{print $1"\t"$1}' | sort -n > $MERGED_BATCH.txt

# (OPTIONAL: REPEAT GENOTYPE CHECK AT MERGED LIBRARY LEVEL)
# SKIPPING FOR NOW
#
# mkdir -p $PILEUP_PATH/$MERGED_BATCH/
# logfile="$PILEUP_PATH/$MERGED_BATCH/run_count_alleles_slurm_%A_%a.log"
# run_batch_task.sh\
#  --array\
#  --batch=$MERGED_BATCH\
#  --sbatchopt="-o $logfile $USENODE"\
#  --task='count_alleles_rna.sh --sample=$sampleID --outpath='$PILEUP_PATH/$MERGED_BATCH/' --header '$STAR_PATH/$MERGED_BATCH/'$sampleID/Aligned.sortedByCoord.out.bam'
#  
# Run the genotyping script on all samples in batch
# sbatch -o genotype_samples.Rout $USENODE genotype_samples.R $MERGED_BATCH


# TO DO: Combine NTR steps into a single one (include incorporation of Wen's NTRs and FlyBase 6 ncRNAs)
# Get rid of the reproducibility criteria, just use full merge
# However, it would be worth counting, for each NTR, how many LINES have overlapping exons from the individual cufflinks transcriptomes


# --- STEP 10: Assemble novel transcriptome --- #

# ----- STEP 10A: Validate Previous NTRs ----- #

# NOTE: This part is only necessary on Baseline, Wen's NTRs will be incorporated into the final NTR set here, which can be reused on other conditions

# Extract just the NTRs from the GTF file Wen created for Tiling Array project:
Rscript Wen_NTR_split_strands.R
NTRGTF="Wen.TilingArray.NTR.both.strands.gtf"

# Run HTSeq-count on this new GTF to determine which models are supported by current data and what strand each transcript is on
for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $HTSEQ_COUNT_PATH/$BATCH
  logfile="$HTSEQ_COUNT_PATH/$BATCH/run_htseq_batch_star_WenNTR_slurm_%A_%a.log"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="$USENODE -o $logfile -c 1" --task='run_htseq.sh --sample=$sampleID --batch=$argBatch --gff='$NTRGTF' --suffix=STAR_WenNTR_counts'
done

# Compile into sample level aggregated count table:
sbatch $USENODE -o build_count_table_WenNTR.Rout ~/Tools/DGRPseq/build_count_table.R STUB=STAR_WenNTR OUTPUT=combined_samples_WenNTR_counts.txt

# Run the final strand assignment script
# and toss out any NTR model here that doesn't have enough reads in the current RNA-seq to support it
Rscript Wen_NTR_assign_strands.R &> Wen_NTR_assign_strands.Rout
# Output file: Wen.TilingArray.NTR.stranded.gtf

# ----- STEP 10B: Assemble individual sample transcriptomes with Cufflinks ----- #

mkdir -p $ASSEMBLED_PATH/$MERGED_BATCH
logfile="$ASSEMBLED_PATH/$MERGED_BATCH/run_cufflinks_star_slurm_%A_%a.log"
./run_batch_task.sh\
 --array\
 --batch=$MERGED_BATCH\
 --sbatchopt="-o $logfile $USENODE -c $DEFAULT_THREAD_COUNT"\
 --task='./run_cufflinks.sh --batch=$argBatch --sample=$sampleID --mode=denovo --max-bundle-frags=10000000000 --max-bundle-length=100000000'

# [CONFIRMED] This should show no difference when complete!
diff <(ls cufflinks/merged_samples/ | grep -v '[.]log' | sort -n) <(cut -f 1 merged_samples.txt) | more

# ----- [DEPRECATED] STEP 10C: Identify reproducible transcripts across biological replicates --- #

# DEPRECATED
# BSAMPNUM=`cut -f 1 ${MERGED_BATCH}.txt | sed 's/[0-9]$//' | sort | uniq | wc -l`

# DEPRECATED
# Identify reproducible transcript regions for each replicate pair
# mkdir -p $ASSEMBLED_PATH/reproducible
# logfile="$ASSEMBLED_PATH/reproducible/extract_reproducible_transcript_%A_%a.log"
# sbatch --array=1-$BSAMPNUM -o $logfile $USENODE extract_reproducible_transcripts_batch.sh $MERGED_BATCH $ASSEMBLED_PATH/reproducible

# NOTE: This new approach (determining reproducibility at the exon level) found ~5k fewer reproducible transcripts overall
# and reduced total novel reproducible transcripts by ~50%
# One possible cause is faulty splice junctions resulting from mismatches at the ends of reads
# These create erroneous splice junctions and additional exons that may appear to be non-reproducible even for transcripts that otherwise are reproducible
# It's possible that cuffmerge can correct this, or at least identify both versions of a transcript model
# So perhaps I *should* run cuffmerge on ALL transcripts, then determine reproducibility based on overlap with exon region intersects from replicate pairs

# Merge across samples and count novel classes of interest
# TO DO: Rename this script now that it just does the merge portion
# TO DO: Include cuffmerge step on transcripts_reproducible.gtf files to get final transcript models
# This script currently just estimates the number of novel features based on merging transcript regions in BED format
# DEPRECATED
# sbatch -o merge_reproducible_transcripts_%j.log $USENODE merge_reproducible_transcripts.sh $MERGED_BATCH $ASSEMBLED_PATH/reproducible

# ----- STEP 10D: Run cuffmerge and cuffcompare to get unified transcriptome ----- #

# Trying this two ways for now:

# I) [DEPRECATED] Run on transcripts_reproducible.gtf files

# mkdir -p $CUFFMERGE_PATH/${MERGED_BATCH}_reproducible
# logfile="$CUFFMERGE_PATH/${MERGED_BATCH}_reproducible/run_cuffmerge_slurm_%j.log"
# sbatch -c 10 -o "$logfile" $USENODE run_cuffmerge.sh --input=$ASSEMBLED_PATH/$MERGED_BATCH --filename=transcripts_reproducible.gtf --output=$CUFFMERGE_PATH/${MERGED_BATCH}_reproducible

# II) Run on original transcripts.gtf files [This is the preferred approach now, no filter for reproducibility]

mkdir -p $CUFFMERGE_PATH/$MERGED_BATCH
logfile="$CUFFMERGE_PATH/$MERGED_BATCH/run_cuffmerge_slurm_%j.log"
sbatch -c 10 -o "$logfile" $USENODE run_cuffmerge.sh --batch=$MERGED_BATCH --prevgtf=Wen.TilingArray.NTR.stranded.gtf

# ----- STEP 10E: Extract novel loci only ----- #

# Count novel gene/transcript categories for each set of cuffmerge/cuffcompare results
# [DEPRECATED] sbatch -o $CUFFMERGE_PATH/${MERGED_BATCH}_reproducible/extract_cuffmerge_novel.log extract_cuffmerge_novel.sh ${MERGED_BATCH}_reproducible
sbatch -o $CUFFMERGE_PATH/$MERGED_BATCH/extract_cuffmerge_novel.log extract_cuffmerge_novel.sh $MERGED_BATCH

# These scripts handle both all and reproducible versions of novel transcriptome

# R script to collect additional stats on novel loci:
./cuffmerge_stats.R $CUFFMERGE_PATH/$MERGED_BATCH

# R script to annotate gene info tables with NTR classes
./annotate_novel_loci.R $CUFFMERGE_PATH/$MERGED_BATCH


# --- STEP 11: Run HTSeq-count on novel features to count STAR-aligned reads --- #

# DEPRECATED - NOT USING REPRODUCIBLE VERSION
# NOVELGTF="$CUFFMERGE_PATH/${MERGED_BATCH}_reproducible/cuffcompare.reproducible.novel.gtf"
#
# for BATCH in ${RUNBATCHES[@]}
# do
#  mkdir -p $HTSEQ_COUNT_PATH/$BATCH
# logfile="$HTSEQ_COUNT_PATH/$BATCH/run_htseq_batch_star_rep_novel_slurm_%A_%a.log"
#  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="$USENODE -o $logfile -c 1" --task='./run_htseq.sh --sample=$sampleID --batch=$argBatch --gff='$NOVELGTF' --suffix=STAR_rep_novel_counts'
# done
#
# # Compile into sample level aggregated count table:
# sbatch $USENODE -o build_count_table_rep_novel.Rout ~/Tools/DGRPseq/build_count_table.R STUB=STAR_rep_novel OUTPUT=combined_samples_rep_novel_counts.txt

# RUN ON THE FULL NOVEL GTF (NO REPRODUCIBLE FILTER)
NOVELGTF="$CUFFMERGE_PATH/${MERGED_BATCH}/cuffcompare.all.novel.gtf"

for BATCH in ${RUNBATCHES[@]}
do
  mkdir -p $HTSEQ_COUNT_PATH/$BATCH
  logfile="$HTSEQ_COUNT_PATH/$BATCH/run_htseq_batch_star_all_novel_slurm_%A_%a.log"
  ./run_batch_task.sh --array --batch=$BATCH --sbatchopt="$USENODE -o $logfile -c 1" --task='run_htseq.sh --sample=$sampleID --batch=$argBatch --gff='$NOVELGTF' --suffix=STAR_all_novel_counts'
done

# Compile into sample level aggregated count table:
sbatch $USENODE -o build_count_table_all_novel.Rout ~/Tools/DGRPseq/build_count_table.R STUB=STAR_all_novel OUTPUT=combined_samples_all_novel_counts.txt


# --- LEADING EDGE --- #
# --- U R HERE - 1/10/17 --- #

