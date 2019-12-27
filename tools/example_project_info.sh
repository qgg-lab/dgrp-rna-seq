#!/bin/bash
#
# Project configuration script - DGRP RNA-seq processing pipeline
#
# Author: Logan J. Everett
# Created: 5/6/16
#
# USAGE: Copy this script to your project directory and rename to "project_info.sh"
#
# This script stores project-wide parameters, and should be loaded each time you log in
# to run commands for that project (see example_pipeline.txt).  This file also acts as a 
# header for other scripts in the pipeline.
#
# In general, you should clone a copy of the DGRPseq toolkit using git, and add that copy
# to your path so that you can run the scripts in all projects.
#
# I recommend creating a separate directory for each project, copy this file into that 
# directory, and then customize the setting for your project.
#
# Typically, the only parameters that need to be changed for a new project are: 
# RAW_FASTQ_PATH=path/to/fastq
# And the list of batch names (for each name "batchX" in this list, there should be a
# file "batchX.txt" that matches the format of batch_example.txt
#
# TO DO: It would be useful to create another script that "initializes" a new project
# by directly asking the user for a few of these parameters and then creating this file
# automatically
#

# PROJECT: [LABEL WHAT PROJECT THIS COPY OF THE FILE IS BEING USED FOR!]


# Path to the raw fastq data to start from
RAW_FASTQ_PATH=/home/hiseq2500/FastqFiles/MY_PROJECT

# Batch names
# Each of these should have a corresponding file batch1.txt, batch2.txt
# See batch_example.txt 
declare -a ALLBATCHES=("batch1" "batch2")

# Path to the top-level directory for this project
PROJECT_HOME='.'

# Multi-threading - use 4 cores per task by default
DEFAULT_THREAD_COUNT=4

# Where to store trimmed fastq files
# This is just processed through cutadapt, but with my new rRNA alignment step they're not really ready for reference genome alignment
# This variable was originally called "PROC_FASTQ_PATH" but was changed "TRIMMED_FASTQ_PATH"
# Although scripts should still keep support for the old name for backwards compatibility
TRIMMED_FASTQ_PATH=$PROJECT_HOME/trimmed

# Genome index to align to (for bowtie2/tophat2)
FLY_GENOME=/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/Bowtie2Index/Dmel_r5.57_FB2014_03

# Genome index for STAR
FLY_GENOME_STAR=/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/StarIndex/

# FASTA file of the genome (for sam file processing)
FLY_GENOME_FASTA=/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta

# Drosophila rRNA sequence database (for BWA)
FLY_RIBOSOMAL_RNA=/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/Ribosomal/flyribo

# Microbiome database (for BWA)
MICROBE_DB=/home/ljeveret/Resources/NCBI/DGRP_Microbiome/DGRP_Microbiome

# Repeat database (for BWA)
REPEAT_DB=/home/ljeveret/Resources/RepBase/dmel.repbase-20.01

# Where to store BWA rRNA alignment results
# (Used to be BWA_PATH before additional microbiome filtering step was added)
BWA_RIBO_PATH=$PROJECT_HOME/bwa_rRNA

# Where to store the final rRNA-filtered fastq files:
RIBO_FILTERED_PATH=$PROJECT_HOME/ribofiltered

# Where to store BWA microbiome DB alignment results
BWA_MICROBE_PATH=$PROJECT_HOME/bwa_microbe

# Where to store microbe-filtered fastq files:
MICROBE_FILTERED_PATH=$PROJECT_HOME/microbefiltered

# Where to store BWA RepBase alignment results
BWA_REPEAT_PATH=$PROJECT_HOME/bwa_repeat

# Where to store repeat-filtered fastq files:
REPEAT_FILTERED_PATH=$PROJECT_HOME/repfiltered

# Where to store STAR alignment results
STAR_PATH=$PROJECT_HOME/star
# TO DO: Should update tophat2 alignment script to make sure it's not over-writing STAR results...

# Where to store Trinity assembly results
TRINITY_PATH=$PROJECT_HOME/trinity

# Where to store cufflinks alignment results by default
ASSEMBLED_PATH=$PROJECT_HOME/cufflinks

# Where to store cuffmerge results for final transcriptome annotation
CUFFMERGE_PATH=$PROJECT_HOME/cuffmerge

# HTSeq-Count defaults parameters
HTSEQ_COUNT_MODE="intersection-nonempty"
HTSEQ_COUNT_STRAND="reverse"
HTSEQ_COUNT_PATH=$PROJECT_HOME/htseq

# GFF File with known transcriptome
FLY_KNOWN_GFF=/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-transcriptome-r5.57.gff
