#!/bin/bash
#
# LJE - 11/10/14
# 
# Script to summarize all fastqc results
#
# Usage:
# ./fastqc_summary.sh FASTQCDIR
#
# Where FASTQCDIR is the name of a top-level directory of fastqc results.
# (Typically fastqc or fastqc_trimmed)
# The expected structure is as follows:
# Must contain a subdirectory for each batch
# Each batch directory must have fully extracted fastqc output (not just the individual html and zip files),
# so each batch directory should contain a subdirectory for each sample matching the pattern (sample)_fastqc/
# and containing files summary.txt and fastqc_data.txt
#
# Currently just want to check for:
# 1) Samples that don't pass basic read quality test
# 2) Make sure all samples have same encoding type
# 3) Check overall %GC 
# 4) Check overall sequence length
# 5) Extract all over-represented sequences and collapse (assemble?) into a unique set
#
# TO DO: This should eventually be replaced with a python script that parses through each section of each fastqc_data.txt file
# 		 and puts the data together into a single chart summarizing all samples
#

# Assign the fastqc primary directory to analyze
FASTQCDIR=$1

# --- STEP 1: Check for key failures in summary.txt files --- #
failout=$FASTQCDIR/key_failures.txt

awk -F'\t' '($2 == "Per base sequence quality" || $2 == "Per sequence quality scores" || $2 == "Per base N content" || $2 == "Sequence Length Distribution")' \
$FASTQCDIR/*/*_fastqc/summary.txt \
| awk '($1 == "FAIL")' \
> $failout

failsz=`wc -l $failout | awk '{print $1}'`

if [[ "$failsz" -le "0" ]]
then
  echo 'No key failures - GOOD!'
else
  echo "Key failures to check out:"
  cat $failout
fi
echo ""

# --- STEP 2: Check encoding types - should all be the same --- #
encout=$FASTQCDIR/encoding_types.txt

awk -F'\t' '($1 == "Encoding")' \
$FASTQCDIR/*/*_fastqc/fastqc_data.txt \
| awk -F'\t' '{print $2}' \
| sort \
| uniq -c \
> $encout

encsz=`wc -l $encout | awk '{print $1}'`

if [[ "$encsz" -eq "1" ]]
then
  echo 'All fastq files are same encoding type - GOOD!'
else
  echo "Some fastq files have alternate encoding types:"
  echo -e "#File Count\t#Encoding Type"
  cat $encsz
fi
echo ""

# --- STEP 3: Check overall %GC --- #
gctmp=$FASTQCDIR/percent_gc.tmp
gcout=$FASTQCDIR/percent_gc.txt

echo -n "" > $gctmp
for qcfile in $FASTQCDIR/*/*_fastqc/fastqc_data.txt
do
  batchname=`echo "$qcfile" | awk -F'/' '{print $2}'`
  sampname=`echo "$qcfile" | awk -F'/' '{print $3}' | sed 's/_fastqc$//g'`
  percgc=`awk -F'\t' '($1 == "%GC")' $qcfile | awk -F'\t' '{print $2}'`
  echo -e "$batchname\t$sampname\t$percgc" >> $gctmp
done

sort -k 3n $gctmp > $gcout
rm $gctmp

echo -e "Summary of GC Content:\n%GC\tSample Count"
awk -F'\t' '{print $3}' $gcout | uniq -c | awk '{print $2"\t"$1}'
echo ""


# --- STEP 4: Check sequence length distribution --- #
seqtmp=$FASTQCDIR/sequence_length.tmp
seqout=$FASTQCDIR/sequence_length.txt

echo -n "" > $seqtmp
for qcfile in $FASTQCDIR/*/*_fastqc/fastqc_data.txt
do
  batchname=`echo "$qcfile" | awk -F'/' '{print $2}'`
  sampname=`echo "$qcfile" | awk -F'/' '{print $3}' | sed 's/_fastqc$//g'`
  seqlen=`awk -F'\t' '($1 == "Sequence length")' $qcfile | awk -F'\t' '{print $2}'`
  echo -e "$batchname\t$sampname\t$seqlen" >> $seqtmp
done

sort -k 3n $seqtmp > $seqout
rm $seqtmp

echo -e "Summary of Sequence Length:\nLength\tSample Count"
awk -F'\t' '{print $3}' $seqout | uniq -c | awk '{print $2"\t"$1}'
echo ""


# --- STEP 5: Over-represented Sequences --- #
knownout=$FASTQCDIR/overrep_known.txt
nohitout=$FASTQCDIR/overrep_nohit.txt

# Extract the over-represented sequences that match a known contaminant
awk -F'\t' '{ if ($1 ~ /^[ACGTN]+$/ && NF == 4 && $4 != "No Hit") print $1"\t"$4 }' $FASTQCDIR/*/*_fastqc/fastqc_data.txt \
| sort \
| uniq \
> $knownout


# Extract the over-represented sequences with no hit
awk -F'\t' '{ if ($1 ~ /^[ACGTN]+$/ && NF == 4 && $4 == "No Hit") print $1"\t"$4 }' $FASTQCDIR/*/*_fastqc/fastqc_data.txt \
| sort \
| uniq \
> $nohitout

knowncnt=`wc -l $knownout | awk '{print $1}'`
nohitcnt=`wc -l $nohitout | awk '{print $1}'`

echo "$knowncnt known sequences are over-represented, see $knownout for details"
echo "$nohitcnt other sequences are over-represented, see $nohitout for details"
echo ""

# TO DO: Find a tool for doing sequence assembly, apply to the No Hit sequences
