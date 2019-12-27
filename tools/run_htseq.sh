#!/bin/bash
# NOTE: Do not set -e here, or will crash when checking samtools version, haven't found a good work around for this yet

# LJE - 11/13/14
#
# run_htseq.sh - Run htseq-count on a single sample
# This version has been stripped of batch functionality - use in combination with run_batch_task.sh
#
# Usage:
# ./run_htseq.sh [Options] [--batch=batchname] --sample=sampleID
#
# NOTE: htseq-count does NOT support multi-threading, so there is no --threads option, ignores DEFAULT_THREAD_COUNT
#
# Options:
#  --sample=sampleID - Refers to the primary directory name where the aligned read BAM file resides, and is also used to name the output file
#		Ultimately, the full input path is: (inputpath)/(batchname)/(sampleID)/(bamfile)
#  --input=path - Path to the input directory, default is $STAR_PATH, if defined in project_info.pl, otherwise defaults to ./
#  --bamfile=string - The standard bam file name in the alignment directory, defaults to "Aligned.sortedByCoord.out.bam"
#  --output=path - Path to the output directory, default is $HTSEQ_COUNT_PATH, as defined in project_info.sh (if missing, defaults to "htseq"
#		If HTSeq-Count output for this sample exists already, this script will quit with an error
#  --suffix=string - Saves the output to (sampleID)_(suffix).txt, default is "counts"
#  --batch=batchname - Convenience option, adds a matching subdirectory to both the input and output paths
#  --gff=path/to/gff - File to use for gene annotations, default is taken from FLY_KNOWN_GFF, if that's also missing this will cause an error
#  --mode=(union|intersection-strict|intersection-nonempty) - Assignment mode for HTSeq-Count, default is intersection-nonempty, or whatever is specified in HTSEQ_COUNT_MODE
#  --strand=(no|yes|reverse) - Strand setting for HTSeq-Count, default is no or whatever is specified in HTSEQ_COUNT_STRAND
#  --overwrite - Over-write any existing output files instead of causing an error (will still warn)
#  --dryrun - Parameter indicating that task should not actually be run, but any pre-task checks should be run,
#		and the final formatted command should be echo'd. This is useful for testing and debugging.
#		Basically, if dryrun is specified, there should be no output or file manipulation (code below must honor this though)
#

# -- COMMAND LINE PARAMETERS -- #

# parse command line options
# NOTE: all variables specified here start with "arg"
args=`getopt -o "s:i:f:o:b:g:m:d" -l "sample:,input:,bamfile:,output:,suffix:,batch:,gff:,mode:,strand:,overwrite,dryrun" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -s|--sample)
      argSample=$2
      shift 2;;

    -i|--input)
      argInput=$2
      shift 2;;
    
    -f|--bamfile)
      argBamFile=$2
      shift 2;;
    
    -o|--output)
      argOutput=$2
      shift 2;;
    
    --suffix)
      argSuffix=$2
      shift 2;;
    
    -b|--batch)
      argBatch=$2
      shift 2;;
    
    -g|--gff)
      argGFF=$2
      shift 2;;
    
    -m|--mode)
      argMode=$2
      shift 2;;
    
    --strand)
      argStrand=$2
      shift 2;;
    
    --overwrite)
      argOverwrite=1
      shift 1;;
    
    -d|--dryrun)
      argDryRun=1
      shift 1;;
    
    --)
      shift
      break;;
  esac
done

# Attempt to load project info
if [[ -e project_info.sh ]]
then
  source project_info.sh
  echo "Loaded project_info.sh"
else
  echo "No project_info.sh found, proceeding without it"
  # TO DO: For many scripts, may want to require this, in which case these two lines should be used instead:
  # echo "ERROR: No project_info.sh found, cannot proceed!"
  # exit 1
fi

# Check each parameter
# Specify required params here and set defaults for others

#  --sample - Required
if [[ ! $argSample ]]
then
  echo "ERROR: Must specify --sample parameter"
  exit 1
fi

#  --input=string - Check for $STAR_PATH if not defined, or default to ./
#			Directory must exist
if [[ ! $argInput ]]
then
  if [[ $STAR_PATH ]] 
  then
    argInput=$STAR_PATH
    echo "Inferring --input=$argInput based on STAR_PATH in project_info.sh"
  else
    argInput="."
    echo "Setting --input=$argInput by default"
  fi
fi

# --bamfile=string - defaults to "Aligned.sortedByCoord.out.bam"
if [[ ! $argBamFile ]]
then
  argBamFile="Aligned.sortedByCoord.out.bam"
  echo "Setting --bamfile=$argBamFile by default"
fi

# --output, if not specified check for HTSEQ_COUNT_PATH, or default to "htseq"
if [[ ! $argOutput ]]
then
  if [[ $HTSEQ_COUNT_PATH ]]
  then
    argOutput=$HTSEQ_COUNT_PATH
    echo "Inferring --output=$argOutput based on HTSEQ_COUNT_PATH in project_info.sh"
  else
    argOutput="htseq"
    echo "Setting --output=$argOutput by default"
  fi
fi

#  --suffix=string - Saves the output to _suffix.txt, default is "counts"
if [[ ! $argSuffix ]]
then
  argSuffix="counts"
  echo "Using default --suffix=$argSuffix"
fi

# --batch - Optional, just append to argInput and argOutput
if [[ $argBatch ]]
then
  argInput=$argInput/$argBatch
  argOutput=$argOutput/$argBatch
  echo "Treating sample as part of batch $argBatch"
  echo "Batch input path: $argInput"
  echo "Batch output path: $argOutput"
fi

# Check that input directory for the current sample exists
if [[ ! -d $argInput/$argSample ]]
then
  echo "ERROR: Input directory does not exist"
  echo "MISSING: $argInput/$argSample"
  exit 1
fi

#  --gff - If not specified, check for FLY_KNOWN_GFF, or throw an error
#		Specified file must exist
if [[ ! $argGFF ]]
then
  if [[ $FLY_KNOWN_GFF ]]
  then
    argGFF=$FLY_KNOWN_GFF
    echo "Inferred --gff=$argGFF from FLY_KNOWN_GFF"
  else
    echo "ERROR: Must specify GFF file with --gff parameter or FLY_KNOWN_GFF environment variable"
    exit 1
  fi
fi

if [[ ! -e $argGFF ]]
then
  echo "ERROR: Missing required file $argGFF"
  exit 1
fi

#  --mode - default to HTSEQ_COUNT_MODE or "union", must be one of (union|intersection-strict|intersection-nonempty)
if [[ ! $argMode ]]
then
  if [[ $HTSEQ_COUNT_MODE ]]
  then
    argMode="$HTSEQ_COUNT_MODE"
    echo "Inferred --mode=$argMode from HTSEQ_COUNT_MODE"
  else
    argMode="union"
    echo "Using default --mode=$argMode"
  fi
fi

if [[ ("$argMode" != "union") && ("$argMode" != "intersection-strict") && ("$argMode" != "intersection-nonempty") ]]
then
  echo "ERROR: --mode=$argMode is not valid, must be one of (union|intersection-strict|intersection-nonempty)"
  exit 1
fi

#  --strand - default to HTSEQ_COUNT_STRAND or no, must be one of (no|yes|reverse)
if [[ ! $argStrand ]]
then
  if [[ $HTSEQ_COUNT_STRAND ]]
  then
    argStrand=$HTSEQ_COUNT_STRAND
    echo "Inferred --strand=$argStrand from HTSEQ_COUNT_STRAND"
  else
    argStrand="no"
    echo "Using default --strand=$argStrand"
  fi
fi

if [[ ($argStrand != "no") && ($argStrand != "yes") && ($argStrand != "reverse") ]]
then
  echo "ERROR: --strand=$argStrand is not valid, must be one of (no|yes|reverse)"
  exit 1
fi

# --dryrun - Just report if this is the case
if [[ $argDryRun ]]
then
  echo "---- DRY RUN ----"
  echo "Tasks will not actually run, there will be no output!"
fi


# TO DO: Check for dependencies and output the versions, e.g.
dependencies=("htseq-count")
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
    else
      $tool --version
    fi
    echo ""
  fi
done


# -- MAIN TASK -- #

# Build and run the main htseq-count command
# Prints the cmd to run, then runs it (stderr from cmd goes to terminal/shell)

# Generate full path to aligned read BAM file
alignedBam=$argInput/$argSample/$argBamFile

# Make sure the file exists
if [[ ! -e $alignedBam ]]
then
  echo "ERROR: Missing $alignedBam"
  exit 1
else
  echo "Counting aligned reads from $alignedBam"
fi
  
# Specify Output file, check if it exists already
countFile=$argOutput/$argSample"_"$argSuffix".txt"
if [[ -e $countFile ]]
then
  if [[ $argOverwrite ]]
  then
    echo "WARNING: HTSeq-count output file already exists:"
    echo "$countFile"
    echo "Will overwrite this file due to user-specified --overwrite parameter"
  else
    echo "ERROR: HTSeq-count output file already exists:"
    echo "$countFile"
    echo "Remove or rename this file before running this task"
    exit 1
  fi
else
  echo "Count data will be in: $countFile"
fi

echo ""
echo "Running: htseq-count -f bam -s $argStrand -m $argMode $alignedBam $argGFF > $countFile"

if [[ ! $argDryRun ]]
then
  mkdir -p $argOutput
  
  echo -ne "Started at:\t"
  date
  echo ""
    
  htseq-count -f bam -s $argStrand -m $argMode $alignedBam $argGFF > $countFile
    
  echo -ne "Finished at:\t"
  date
  echo ""
fi

echo "Entire script completed successfully."
