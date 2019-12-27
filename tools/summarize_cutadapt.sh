#!/bin/bash
#
# LJE - 11/11/14
#
# Pulls key information out of the log file from a single cutadapt run
# Outputs the results to a single tab-delimited line, so if run on a batch of files it will form a tab-delimited table
#
# Usage:
# ./summarize_cutadapt.sh file1.log [file2.log ...]
#
# Output format:
# fastqfile, batch, totalReads, shortReads, emptyReads
#
# Where fastqfile, batch refer to the particular sample
# totalReads is the "Processed Reads", or total prior to any trimming
# shortReads is how many were filtered for being too short
# emptyReads is how many where the full length (125) was removed

echo -e "#FastqFile\tBatch\tTotalReads\tShortReads\tEmptyReads"

for logfile in $@
do
  fastqFile=`grep '^Running:' $logfile | awk '{print $8}' | awk -F'/' '{print $4}'`
  batchName=`grep '^Running:' $logfile | awk '{print $8}' | awk -F'/' '{print $3}'`
  totalReads=`grep 'Processed reads:' $logfile | grep -o '[0-9]*'`
  shortReads=`grep 'Too short reads:' $logfile | awk '{print $4}'`
  emptyReads=`grep '^[1-9]' $logfile | tail -n 1 | awk '{print $2}'`
  echo -e "$fastqFile\t$batchName\t$totalReads\t$shortReads\t$emptyReads"
done
