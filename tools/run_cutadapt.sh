#!/bin/bash
set -e

# LJE - 3/24/15
#
# run_cutadapt.sh - Script for running cutadapt on sequence files for DGRP Baseline RNA-seq project
# This no longer handles whole batches or slurm job submission, use in combination with run_batch_task.sh for that
#
# Usage:
#  ./run_cutadapt.sh [options] --sample=sampleID
#
# OPTIONS:
#  --sample=sampleID (Required) - refers to the input file name, except for the suffix (also excludes the _1, _2 for paired sequence files)
#	e.g. --sample=mydata will look for mydata.fastq, and if --paired is also specified then it will look for mydata_1.fastq and mydata_2.fastq
#	This script will auto-detect if the expected input files have been compressed.
#  --batch=batchName - Specifies batch subdirectory for the output
#  --input=inputpath - Specifies the path to the input file(s). Defaults to RAW_FASTQ_PATH, specified in project_info.sh
#  --output=outpath - Specifies the path to store the output file(s). 
#	Defaults to $[TRIMMED|PROC]_FASTQ_PATH/$BATCH if these are set as environment variables, or ./trimmed if nothing else specified
#	Output files are stored in this path as $sampleID(_1/2)_trimmed.fastq.gz
#	NOTE: Script will quit with error if the output file exists already!
#  --paired - Tells the script to process paired end sequence data (processes single end by default)
#  --dryrun - Just test that required input files are present and output files are not present, don't actually run anything
#
# A few fixed parameter adjustments have been made here:
# 1) Getting rid of quality filtering - tophat2 and STAR handle that themselves
# 2) Increasing min-length to 50 - shouldn't attempt to align things shorter than that as it will probably inflate ambiguous mapping %
#

# parse command line options
args=`getopt -o "s:b:i:o:p" -l "sample:,batch:,input:,output:,paired,dryrun" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -s|--sample)
      argSample=$2
      shift 2;;
    
    -b|--batch)
      argBatch=$2
      shift 2;;
    
    -i|--input)
      argInput=$2
      shift 2;;
    
    -o|--output)
      argOutput=$2
      shift 2;;
    
    -p|--paired)
      argPaired=1
      shift 1;;
    
    --dryrun)
      argDryRun=1
      shift 1;;
    
    --)
      shift
      break;;
  esac
done

# Attempt to load project info for certain defaults
if [[ -e project_info.sh ]]
then
  source project_info.sh
  echo "Loaded project_info.sh"
else
  echo "No project_info.sh, attempting to proceed without it"
fi

# --- Options --- #

if [[ ! $argSample ]]
then
  echo "ERROR: Must specify a sample ID!"
  exit 1
fi

# If --input isn't set, check RAW_FASTQ_PATH
if [[ ! $argInput ]]
then
  if [[ $RAW_FASTQ_PATH ]]
  then
    argInput=$RAW_FASTQ_PATH
    echo "Setting --input=$argInput based on RAW_FASTQ_PATH in project_info.sh."
  else
    echo "ERROR: Must specify input path using --input parameter or set RAW_FASTQ_PATH prior to calling this script!"
    exit 1
  fi
fi

# Make sure argInput is a valid directory
if [[ ! -d $argInput ]]
then
  echo "ERROR: $argInput is not a valid directory!"
  exit 1
fi

# If --output isn't set, check TRIMMED_FASTQ_PATH, then check PROC_FASTQ_PATH (older project_info.sh files)
if [[ ! $argOutput ]]
then
  if [[ $TRIMMED_FASTQ_PATH ]]
  then
    argOutput=$TRIMMED_FASTQ_PATH
    echo "Inferring --output=$argOutput based on TRIMMED_FASTQ_PATH in project_info.sh"
  elif [[ $PROC_FASTQ_PATH ]]
  then
    argOutput=$PROC_FASTQ_PATH
    echo "Inferring --output=$argOutput based on PROC_FASTQ_PATH in project_info.sh"
  else
    argOutput="trimmed"
    echo "Inferring --output=$argOutput by default"
  fi
  # Now check for $argBatch or $BATCH variables - append if either is present
  if [[ $argBatch ]]
  then
    argOutput="$argOutput/$argBatch"
  fi
fi

# Make sure output directory exists
mkdir -p $argOutput

# -- MAIN PROCESSING -- #

# The full paths to the input and output file(s are inferred from variables above
# Prints the cmd to run, then runs it (output goes to stdout/stderr)

if [[ ! $argPaired ]]
then
  fastqin=$argInput/$argSample.fastq
  fastqout=$argOutput/${argSample}_trimmed.fastq.gz
  
  if [[ ! -e $fastqin ]]
  then
    if [[ -e $fastqin.gz ]]
    then
  	  fastqin=$fastqin.gz
    else
  	  echo "ERROR: Missing $fastqin(.gz)"
  	  exit 1
    fi
  fi
  
  if [[ -e $fastqout || -e $fastqout.gz ]]
  then
    echo "ERROR: $fastqout(.gz) already exists, remove before re-running this script."
    exit 1
  fi
  
  echo "Running: cutadapt -a TruSeq=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 50 -o $fastqout $fastqin"
  
  if [[ ! $argDryRun ]]
  then
    cutadapt -a TruSeq=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 50 -o $fastqout $fastqin
  fi
else
  # PAIRED SEQUENCE FILES
  echo "Processing paired sequence files."
  echo ""
  
  # Note: for some reason outputting directly to gzipped files doesn't work with paired files
  
  # UPDATE - 9/7/16
  # Make compatible with newer fastq file formats
  if [[ $argSample =~ [_]R1[_] ]]
  then
    leftstub=${argSample}
    # NOTE: Reversing before sed guarantees that the LAST instance of _R1_ gets replace with _R2_ (in case _R1_ also part of sample name)
    rightstub=`echo $argSample | rev | sed 's/[_]1R[_]/_2R_/' | rev`
  else
    leftstub=$argSample"_1"
    rightstub=$argSample"_2"
  fi
  
  leftin=$argInput/${leftstub}.fastq
  rightin=$argInput/${rightstub}.fastq
  lefttmp=$argOutput/${leftstub}_cutadapt.tmp.fastq
  righttmp=$argOutput/${rightstub}_cutadapt.tmp.fastq
  leftout=$argOutput/${leftstub}_trimmed.fastq
  rightout=$argOutput/${rightstub}_trimmed.fastq
  
  if [[ ! -e $leftin ]]
  then
    if [[ -e $leftin.gz ]]
    then
  	  leftin=$leftin.gz
    else
  	  echo "ERROR: Missing $leftin(.gz)"
  	  exit 1
    fi
  fi
  if [[ ! -e $rightin ]]
  then
    if [[ -e $rightin.gz ]]
    then
  	  rightin=$rightin.gz
    else
  	  echo "ERROR: Missing $rightin(.gz)"
  	  exit 1
    fi
  fi
  if [[ -e $lefttmp || -e $righttmp ]]
  then
    echo "ERROR: $lefttmp and/or $righttmp already exists, remove before running this script."
    exit 1
  fi
  if [[ -e $leftout || -e $leftout.gz || -e $rightout || -e $rightout.gz ]]
  then
    echo "ERROR: $leftout(.gz) and/or $right.out(.gz) already exists, remove before running this script."
    exit 1
  fi
  
  echo "Running: cutadapt -a TruSeq=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 50 -o $lefttmp -p $righttmp $leftin $rightin"
  
  if [[ ! $argDryRun ]]
  then
    cutadapt -a TruSeq=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 50 -o $lefttmp -p $righttmp $leftin $rightin
  fi
  
  echo ""
  echo "Running: cutadapt -a TruSeqRev=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 50 -o $rightout -p $leftout $righttmp $lefttmp"
  
  if [[ ! $argDryRun ]]
  then
    cutadapt -a TruSeqRev=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 50 -o $rightout -p $leftout $righttmp $lefttmp
    
    echo ""
    echo "Temp file stats:"
    wc -l $lefttmp
    wc -l $righttmp
    wc -l $leftout
    wc -l $rightout
    echo "Removing $lefttmp and $righttmp"
    rm $lefttmp $righttmp
    echo "Compressing $leftout and $rightout"
    gzip $leftout $rightout
  fi
fi
