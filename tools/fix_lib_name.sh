#!/bin/bash
#
# LJE - 9/1/16
#
# Script to fix library names for samples that were mixed up
# It is best to copy this script into your project directory and tailor to the current state of the project at the time when you are correcting sample names
# Also a good idea to make another script (e.g. fix_all_libs.sh) that runs the script for ALL sample name corrections
# Running the script again for samples that have already been corrected should not cause any problems
#
# IMPORTANT: You must correct the original fastq file (in RAW_FASTQ_PATH) BEFORE running this script!
# You should also change the corresponding line in your batch file first (this script will verify that you did it right)
#
# USAGE:
# ./fix_lib_name.sh BATCH CURNAME NEWNAME
#
# You should run this before cleaning up any intermediate files!
# If the script finds any other files still using old sample name, clean them up manually, or add them to the renameFile commands below
# AFTER correcting all mixed up file names, you need to rerun the steps that build fastq.gz.stats files, 
# then QC tables, sample info tables, etc. to propagate the corrected sample names
#
# Once all file names have been corrected, you should re-run genotype_samples.R on all batches/samples and make sure the relevant errors go away
#

BATCH=$1
CURNAME=$2
NEWNAME=$3

if [[ "$CURNAME" == "$NEWNAME" ]]
then
  echo "ERROR: NEWNAME must be different from CURNAME"
  exit 1
else
  echo "Renaming library $CURNAME to $NEWNAME"
fi

# LOAD PROJECT-WIDE PARAMS
./project_info.sh

# FIRST: check that original library has been renamed, if not quit
CURSRC="$RAW_FASTQ_PATH/$CURNAME.fastq.gz"
NEWSRC="$RAW_FASTQ_PATH/$NEWNAME.fastq.gz"

if [[ -e $CURSRC ]]
then
  if [[ -e $NEWSRC ]]
  then
    echo "ERROR: $CURSRC and $NEWSRC both exist"
    ls -lh $CURSRC $NEWSRC
    md5sum $CURSRC $NEWSRC
    exit 1
  else
    echo "ERROR: Rename $CURSRC to $NEWSRC first"
    echo "Run this command with necessary permissions:"
    echo "mv $CURSRC $NEWSRC"
    exit 1
  fi
else
  if [[ -e $NEWSRC ]]
  then
    echo "OK: $CURSRC already renamed to $NEWSRC"
  else
    echo "ERROR: Neither $CURSRC nor $NEWSRC exists, mistake in library names?"
    exit 1
  fi
fi

# --- Sub-routine to rename a file if it hasn't been renamed already --- #

# NOTE: This is safe for files OR directories
# For directories, make sure not to have a trailing / on either cur or new file
function renameFile {
  CURFILE=$1
  NEWFILE=$2
  if [[ -e $CURFILE ]]
  then
    if [[ -e $NEWFILE ]]
    then
      echo "ERROR: $CURFILE and $NEWFILE both exist"
      ls -lh $CURFILE $NEWFILE
      md5sum $CURFILE $NEWFILE
      exit 1
    else
      echo "$CURFILE -> $NEWFILE"
      mv $CURFILE $NEWFILE
    fi
  else
    if [[ -e $NEWFILE ]]
    then
      echo "$CURFILE already moved to $NEWFILE"
    else
      echo "WARNING: Could not find $CURFILE or $NEWFILE"
    fi
  fi
}

# If any files are no longer used in project, you can comment out the relevant lines below to avoid errors

# FASTQC RESULTS:
renameFile fastqc/$BATCH/${CURNAME}_fastqc.zip fastqc/$BATCH/${NEWNAME}_fastqc.zip
renameFile fastqc/$BATCH/${CURNAME}_fastqc fastqc/$BATCH/${NEWNAME}_fastqc
renameFile fastqc_repfiltered/$BATCH/${CURNAME}_filtered_fastqc.zip fastqc_repfiltered/$BATCH/${NEWNAME}_filtered_fastqc.zip
renameFile fastqc_repfiltered/$BATCH/${CURNAME}_filtered_fastqc fastqc_repfiltered/$BATCH/${NEWNAME}_filtered_fastqc

# BWA RESULTS (Directories AND fastq.gz files):
renameFile $BWA_RIBO_PATH/$BATCH/${CURNAME} $BWA_RIBO_PATH/$BATCH/${NEWNAME}
# The directory is renamed, but the fastq.gz file is not:
renameFile $BWA_RIBO_PATH/$BATCH/${NEWNAME}/${CURNAME}_mapped.fastq.gz $BWA_RIBO_PATH/$BATCH/${NEWNAME}/${NEWNAME}_mapped.fastq.gz
renameFile $BWA_MICROBE_PATH/$BATCH/${CURNAME} $BWA_MICROBE_PATH/$BATCH/${NEWNAME}
# The directory is renamed, but the fastq.gz file is not:
renameFile $BWA_MICROBE_PATH/$BATCH/${NEWNAME}/${CURNAME}_mapped.fastq.gz $BWA_MICROBE_PATH/$BATCH/${NEWNAME}/${NEWNAME}_mapped.fastq.gz
renameFile $BWA_REPEAT_PATH/$BATCH/${CURNAME} $BWA_REPEAT_PATH/$BATCH/${NEWNAME}
# The directory is renamed, but the fastq.gz file is not:
renameFile $BWA_REPEAT_PATH/$BATCH/${NEWNAME}/${CURNAME}_mapped.fastq.gz $BWA_REPEAT_PATH/$BATCH/${NEWNAME}/${NEWNAME}_mapped.fastq.gz

# Filtered fastq.gz files:
renameFile $TRIMMED_FASTQ_PATH/$BATCH/${CURNAME}_trimmed.fastq.gz $TRIMMED_FASTQ_PATH/$BATCH/${NEWNAME}_trimmed.fastq.gz
renameFile $RIBO_FILTERED_PATH/$BATCH/${CURNAME}_filtered.fastq.gz $RIBO_FILTERED_PATH/$BATCH/${NEWNAME}_filtered.fastq.gz
renameFile $MICROBE_FILTERED_PATH/$BATCH/${CURNAME}_filtered.fastq.gz $MICROBE_FILTERED_PATH/$BATCH/${NEWNAME}_filtered.fastq.gz 
renameFile $REPEAT_FILTERED_PATH/$BATCH/${CURNAME}_filtered.fastq.gz $REPEAT_FILTERED_PATH/$BATCH/${NEWNAME}_filtered.fastq.gz

# STAR OUTPUT:
renameFile $STAR_PATH/$BATCH/${CURNAME} $STAR_PATH/$BATCH/${NEWNAME}

# Allele Counts
renameFile $PILEUP_PATH/$BATCH/${CURNAME}_counts.txt $PILEUP_PATH/$BATCH/${NEWNAME}_counts.txt

# HTSeq Results: (If correcting before running HTSeq, comment this line out)
renameFile $HTSEQ_COUNT_PATH/$BATCH/${CURNAME}_STAR_counts.txt $HTSEQ_COUNT_PATH/$BATCH/${NEWNAME}_STAR_counts.txt

echo ""

# If any additional files are found, correct them manually, or add the pattern to the sequence of renameFile commands above
echo "Checking for any other files matching $CURNAME ..."
find . -iname '*'$CURNAME'*'
echo ""

# Don't forget to update the batch file!
echo "Checking for old batch file line:"
grep $CURNAME $BATCH.txt
echo "Checking for new batch file line:"
grep $CURNAME $BATCH.txt
echo ""

echo "Finished."
echo ""
