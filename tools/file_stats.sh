#!/bin/bash
set -e

#
# LJE - 2/15/16
#
# Collect basic stats on large files (number of lines, size in bytes, md5sum)
# Typically used for collecting stats on read files like fastq and bam
# 
# Usage:
# ./file_stats.sh data/*.fastq > data/fastq.stats
#

# Loop over parameters - each one should be a file or file pattern
# The inner loop guarantees that file patterns will be expanded
for bigfile in $@
do
  # Collect each stat
  # Read or line count - this depends on file type!
  if [[ "$bigfile" =~ [.]gz$ ]]
  then
    # $bigfile = GZIPPED FILE
    # Pipe through gunzip -c to decompress
    LCOUNT=`gunzip -c $bigfile | wc -l`
  elif [[ "$bigfile" =~ [.][bs]am$ ]]
  then
    # $bigfile = BAM/SAM
    # Use samtools to get read count
    LCOUNT=`samtools view $bigfile | wc -l`
  else
    # Assume $bigfile = PLAIN TEXT (including fastq, etc.)
    # Just count lines
    LCOUNT=`wc -l $bigfile | awk '{print $1}'`
  fi
  # Get file size in bytes
  FSIZE=`du -b $bigfile | awk '{print $1}'`
  # Get MD5 Checksum
  CHECKSUM=`md5sum $bigfile | awk '{print $1}'`
  # Write output line
  echo -e "$bigfile\t$LCOUNT\t$FSIZE\t$CHECKSUM"
done
