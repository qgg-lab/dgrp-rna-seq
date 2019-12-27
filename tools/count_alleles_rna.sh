#!/bin/bash
set -e

# LJE - 8/16/16
#
# Run mpileup and post-processing scripts on one or more bam files to count allele frequencies
# at a pre-determined set of SNPs and flag problematic instances
# NOTE: This script does NOT handle batch/array logic, it is designed to be run as the 
# underlying task for run_batch_task.sh in that context.
# However, it does check project_info.sh for necessary script defaults
# NOTE: Adapted from count_alleles.sh script in Lab_Evolution project
# This version is intended for RNA-seq data, although the logic might not be much different than for DNA
# Can probably combine the scripts in this repository later.
#
# Usage:
# ./count_alleles.sh [Options] file1.bam [file2.bam ...]
#
# OUTPUT: stdout of the mpileup/post-proc chain goes to: outpath/outfile (see params below)
#		  stderr still goes to terminal (or log file from sbatch task)
# 		  If output/outfile exists already, this script will exit with error
# 
# GUIDE TO PARAMETERS:
#  -s | --sample=NAME = Name of sample, only used if --outfile is not specified
#  --input=path - Path to the input directory, default is $STAR_PATH, if defined in project_info.pl, otherwise defaults to ./
#  -f | --outfile=NAME = Name of output file (processed mpileup results)
#	Defaults to (sample)_counts.txt if --sample specified, allele_counts.txt otherwise
#  -p | --outpath=PATH = Path for where to put output file
#	If not specified, defaults to PILEUP_PATH if defined, or . otherwise
#   If outpath does not exist, it will be created before running the main task
#  -t | --trackfile=NAME = Name of track file (I'm not actually sure if anything useful goes in here, it comes out of pileupCov.pl)
#	Default: (outpath)/(sample)_track.out or track_dump.out if not specified
#   Unlike outfile, you will only get a warning if overwriting an existing file here
#  -h | --header = Add a header line to the output
#  -d | --dryrun = Path 
#
# TO DO: Parameterize other options currently pulled from project_info.sh
# TO DO: Option to add a header line
# TO DO: Option to assign shorter name to each input file
#
# REQUIRES: samtools, pileupCov.pl, and filterVariants.pl must be in PATH
#

# get command line options
args=`getopt -o "s:f:p:t:hd" -l "sample:,outfile:,outpath:,trackfile:,header,dryrun" -- "$@"`
echo "Running with parameters: $args"
eval set -- "$args"

# parse arguments
while true;
do
  case $1 in
    -s|--sample)
      argSample=$2
      shift 2;;

	-f|--outfile)
	  argOutFile=$2
	  shift 2;;
      
    -p|--outpath)
      argOutPath=$2
      shift 2;;
    
    -t|--trackfile)
      argTrackFile=$2
      shift 2;;
    
    -h|--header)
      argHeader=1
      shift 1;;
    
    -d|--dryrun)
      argDryRun=1
      shift 1;;
    
    --)
      shift
      break;;
  esac
done

# Remaining parameters are the bam files to process with mpileup
if [[ "${#@}" -eq "0" ]]
then
  echo "ERROR: Must specify one or more input bam files for mpileup"
  exit 1
else
  argInputs=("$@")
fi

# Attempt to load project info first, useful for defaults
if [[ -e project_info.sh ]]
then
  source project_info.sh
  echo "Loaded project_info.sh"
else
  # echo "No project_info.sh found, proceeding without it"
  # TO DO: If all other fields from project_info.sh can be set on command line, then it would be safe to proceed
  echo "ERROR: No project_info.sh found, cannot proceed!"
  exit 1
fi

# Check parameter validity

# --outfile - Defaults to (sample)_counts.txt if --sample specified, allele_counts.txt otherwise
if [[ ! $argOutFile ]]
then
  if [[ $argSample ]]
  then
    argOutFile="${argSample}_counts.txt"
  else
    argOutFile="allele_counts.txt"
  fi
fi

# --outpath defaults to PILEUP_PATH if defined, or . otherwise
if [[ ! $argOutPath ]]
then
  if [[ $PILEUP_PATH ]]
  then
    argOutPath=$PILEUP_PATH
    echo "Inferring --outpath=$argOutPath based on PILEUP_PATH from project_info.sh"
  else
    argOutPath="."
  fi
fi

# --trackfile defaults to (sample)_track.out or track_dump.out
if [[ ! $argTrackFile ]]
then
  if [[ $argSample ]]
  then
    argTrackFile="${argSample}_track.out"
  else
    argTrackFile="track_dump.out"
  fi
fi

# --dryrun - Just alert if set
if [[ $argDryRun ]]
then
  echo ""
  echo "---- DRY RUN ----"
  echo "Tasks will not actually run, there will be no disk output!"
  echo "---- DRY RUN -----"
  echo ""
fi

# Additional parameters (pulled from project_info.pl or go to defaults)
# TO DO: Allow these to be set from the command line as well

# This is the SNP file to do final allele counts for
# TO DO: Confirm this is actually the appropriate file to use for bristle selection project
# Was bristle selection project started from FlyLand population?
if [[ ! $SNP_MAP ]]
then
  SNP_MAP="/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/freeze2.all.snp.bed"
  echo "Using default SNP_MAP = $SNP_MAP"
else
  echo "SNP_MAP = $SNP_MAP"
fi

if [[ ! $FLY_GENOME_FASTA ]]
then
  FLY_GENOME_FASTA="/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta"
  echo "Using default FLY_GENOME_FASTA = $FLY_GENOME_FASTA"
else
  echo "FLY_GENOME_FASTA = $FLY_GENOME_FASTA"
fi

# Check for dependencies and output the versions, e.g.
dependencies=("samtools" "pileupCov.pl" "filterVariants.pl")
echo ""
echo "Checking for required tools: ${dependencies[@]}"
echo ""

for tool in ${dependencies[@]}
do
  toolpath=`which $tool`
  if [[ $toolpath == "" ]]
  then
    echo "ERROR: Missing $tool"
    exit 1
  else
    echo "$toolpath"
    if [[ $tool == "samtools" || $tool == "bcftools" ]]
    then
      $tool 2>&1 | grep '^Version:'
    elif [[ $tool == "htseq-count" ]]
    then
      $tool 2>&1 | tail -n 3
    elif [[ $tool == "pileupCov.pl" || $tool == "filterVariants.pl" ]]
    then
      ls -l $toolpath
    else
      $tool --version
    fi
    if [[ $tool == "STAR" ]]
    then
      echo ""
    fi
    echo ""
  fi
done


# Check if output file exists already (Quit if it does)
outfile="$argOutPath/$argOutFile"
if [[ -e $outfile ]]
then
  echo "ERROR: $outfile exists already"
  echo "Remove or rename this file before running this command"
  exit 1
else
  echo "Saving primary output to: $outfile"
  echo ""
fi

# Also check on the track file - just warn if it exists
if [[ -e $argOutPath/$argTrackFile ]]
then
  echo "WARNING: $argOutPath/$argTrackFile will be overwritten!"
  echo ""
# else
#  echo "Track output: $argOutPath/$argTrackFile"
#  echo ""
fi
# NOTE: NO TRACK OUTPUT WITH -mono 0 PARAM ON pileupCov.pl


# MAIN TASK: Create pileup format table

# TO DO: Add parameters to adjust other options on each of these tools?

echo "RUNNING: samtools mpileup -Q 1 -q 13 -d 10000 -l $SNP_MAP ${argInputs[@]} | pileupCov.pl -minQ 13 -track $argOutPath/$argTrackFile -ref $FLY_GENOME_FASTA -mono 0 | filterVariants.pl -minGB 10 -minTC 10 >> $outfile"

if [[ ! $argDryRun ]]
then

  mkdir -p $argOutPath

  echo ""
  echo -ne "MPILEUP START:\t"
  date
  echo ""

  # Skipping this filtering step because it doesn't generalize to arbitrary number of samples
  # And also because I want to see when certain SNPs get an error flag across the board
  # awk '($6 == "AF" || $6 == "PASS") && ($8 == "AF" || $8 == "PASS")'
  
  if [[ $argHeader ]]
  then
    echo -ne "CHR\tPOS\tREF\tMAJOR\tMINOR" > $outfile
    for inputFile in "${argInputs[@]}"
    do
      echo -ne "\t${inputFile}.FILTER\t${inputFile}.COUNTS" >> $outfile
    done
    echo "" >> $outfile
  else
    echo -n "" > $outfile
  fi
  
  samtools mpileup -Q 1 -q 13 -d 10000 -l $SNP_MAP ${argInputs[@]} \
  | pileupCov.pl -minQ 13 -track $argOutPath/$argTrackFile -ref $FLY_GENOME_FASTA -mono 0 \
  | filterVariants.pl -minGB 10 -minTC 10 \
  >> $outfile
  
  mpileup_exit_code=$?
  
  echo ""
  echo -ne "MPILEUP FINISH:\t"
  date
  echo ""
  
  if [[ "$mpileup_exit_code" -ne "0" ]]
  then
    echo "ERROR: mpileup command returned non-zero exit code = $mpileup_exit_code"
    exit $mpileup_exit_code
  fi
fi
  



