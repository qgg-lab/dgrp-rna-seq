#!/bin/bash
# LJE - 12/8/15
#
# Loop over all processed fastq files to get line counts 
# Assumes batch subdirectories, and that all fastq files are gzipped
#
# Usage:
# fastq_sizes.sh path
#
# Loops over each file matching path/*/*.fastq.gz, computes line count,
# and outputs two-column format (file <tab> line count) to path/fastq.gz.wcl

mypath=$1
outfile=$mypath/fastq.gz.wcl

echo -n "" > $outfile

for fqgz in $mypath/*/*.fastq.gz
do
   lc=`zcat $fqgz | wc -l`
   echo -e "$fqgz\t$lc" >> $outfile
done

echo "Fastq size calculations complete."
