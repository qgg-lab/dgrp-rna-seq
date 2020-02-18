#!/bin/bash

# Job name:
#SBATCH --job-name=dm-cell-line
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Memory per node:
#SBATCH --mem=8G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=2:00:00

echo "SLURM_NTASKS: $SLURM_NTASKS"
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"

cd /mnt/home/huangw53/qgg/dgrp-cell-line-exp/

~/qgg/software/kallisto-v0.46.1/kallisto quant -i all.idx -o quant/"$srr".quant /mnt/home/huangw53/scratch/dmCellLines/fastq/"$srr"_1.fastq /mnt/home/huangw53/scratch/dmCellLines/fastq/"$srr"_2.fastq > quant/"$srr".quant.log 2>&1

sacct -j $SLURM_JOB_ID --units=G -P --format="Account%15,JobID,JobName%15,CPUTime,Elapsed,MaxRSS,Start,End,NodeList,NTasks,Partition%20,NCPUS,state"
