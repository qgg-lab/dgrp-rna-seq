#!/bin/bash
#
# Takes a single parameter which is the directory to run on
# Target dir should contain a file called Aligned.sortedByCoord.out.bam
# Outputs a file called reads_aligned_Dmel_star.txt in target dir
#

BAMFILE=$1/Aligned.sortedByCoord.out.bam
READFILE=$1/reads_aligned_Dmel_star.txt

echo "Extracting read IDs from $BAMFILE to $READFILE"

samtools view -F 4 $BAMFILE | awk -F'\t' '($3 !~ "^U") {print $1}' | sort | uniq > $READFILE
