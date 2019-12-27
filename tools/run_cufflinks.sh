#!/bin/bash
# NOTE: Do not set -e here, or will crash when checking samtools version, haven't found a good work around for this yet

# LJE - 4/2/15
#
# run_cufflinks.sh - Run cufflinks on a single sample
# This version has been stripped of batch functionality - use in combination with run_batch_task.sh
#
# Usage:
# ./run_cufflinks.sh [options] --batch=batchname --sample=sampleID
#
# Options:
#  --sample=sampleID - Refers to the primary directory name where the aligned read BAM file resides, and is also used to name the output directory
#		Ultimately, the full input path is: (inputpath)/(batchname)/(sampleID)/(bamfile)
#  --input=path - Path to the input directory, default is $STAR_PATH, if defined in project_info.pl, otherwise defaults to ./
#  --bamfile=string - The standard bam file name in the alignment directory, defaults to "Aligned.sortedByCoord.out.bam"
#  --output=path - Path to the output directory, default is $ASSEMBLED_PATH, as defined in project_info.sh (if missing, defaults to "cufflinks"
#		The actual output goes in a subdirectory named sampleID in the path indicated here
#		If STAR output for this sample exists already, this script will quit with an error
#  --batch=batchname - Convenience option, adds a matching subdirectory to both the input and output paths
#  --gff=path/to/gff - File to use for gene annotations, default is taken from FLY_KNOWN_GFF, 
#			if that's also missing this will cause an error (unless --mode=denovo)
#  --mode=guide|known|denovo - Controls the -g/-G gff parameter passed to cufflinks
#			guide (Default) - Uses the -g gff option to specify a known transcriptome used to guide the assembly (RABT mode)
#  			known - Uses the -G gff option to specify a known transcriptome to be quantified (no novel transcript discovery)
#			denovo - Suppresses use of a known transcriptome, does assembly de novo
#  --library-type - Passed directly to cufflinks, default is fr-firststrand
#  --max-bundle-frags - Passed directly to cufflinks, default is 1000000 but can be raised if too much stuff in skipped.gtf
#			If not specified, this parameter will not get passed to cufflinks, so it will default to whatever cufflink's default is
#			10000000 is recommended for assembling the ribo_only portion
#  --max-bundle-length - Passed directly to cufflinks, can be raised it too much stuff in skipped.gtf
#			If not specified, this parameter will not get passed to cufflinks, so it will default to whatever cufflink's default is
#  --extraparams - Pass additional parameters to cufflinks that are aren't currently part of the options here
#			Just make sure to protect with single quotes
#  --threads=N - Specifies the number of threads to use, if unspecified, will check SLURM environment, then DEFAULT_THREAD_COUNT, otherwise defaults to 1
#  --dryrun - Parameter indicating that task should not actually be run, but any pre-task checks should be run,
#		and the final formatted command should be echo'd. This is useful for testing and debugging.
#		Basically, if dryrun is specified, there should be no output or file manipulation (code below must honor this though)
#
# UPDATE 2/6/15 - Running into issues where huge portions of each chromosome is getting lumped into one bundle and ends up in skipped.gtf  
# After a bit of investigation, I don't think this is some mega-transcript spanning lots of densely packed genes
# I've tried increasing the --max-bundle-frags and --max-bundle-length parameters quite high
# TODO try adjusting the --max-mle-iterations parameter (default=5000)
# NOTE: with -G option this problem doesn't happen, and things run MUCH faster (4 minutes vs 20 hours with 10 threads)
# So this is clearly a problem when doing novel transcript discovery
# RUNNING in "de novo" mode - no GTF guide at all - also works just fine, so that may be the only way to deal with this...
#

# -- COMMAND LINE PARAMETERS -- #

# parse command line options
# NOTE: all variables specified here start with "arg"
args=`getopt -o "s:i:f:o:b:g:m:t:d" -l "sample:,input:,bamfile:,output:,batch:,gff:,mode:,library-type:,max-bundle-frags:,max-bundle-length:,extraparams:,threads:,dryrun" -- "$@"`
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
        
    -b|--batch)
      argBatch=$2
      shift 2;;
    
    -g|--gff)
      argGFF=$2
      shift 2;;
    
    -m|--mode)
      argMode=$2
      shift 2;;
    
    --library-type)
      argLibraryType=$2
      shift 2;;
    
    --max-bundle-frags)
      argMaxBundleFrags=$2
      shift 2;;
    
    --max-bundle-length)
      argMaxBundleLength=$2
      shift 2;;
    
    --extraparams)
      argExtraParams=$2
      shift 2;;
    
    -t|--threads)
      argThreads=$2
      shift 2;;
        
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

# --sample - Required
if [[ ! $argSample ]]
then
  echo "ERROR: Must specify --sample parameter"
  exit 1
fi

# --input=string - Check for $STAR_PATH if not defined, or default to ./
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

# --output, if not specified check for ASSEMBLED_PATH, or default to "cufflinks"
if [[ ! $argOutput ]]
then
  if [[ $ASSEMBLED_PATH ]]
  then
    argOutput=$ASSEMBLED_PATH
    echo "Inferring --output=$argOutput based on ASSEMBLED_PATH in project_info.sh"
  else
    argOutput=cufflinks
    echo "Setting --output=$argOutput by default"
  fi
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

# --threads - If not specified, check for SLURM environment variables, then DEFAULT_THREAD_COUNT, or default to 1
if [[ ! $argThreads ]]
then
  if [[ $SLURM_JOB_CPUS_PER_NODE ]] 
  then
    argThreads=$SLURM_JOB_CPUS_PER_NODE
    echo "Setting --threads=$argThreads based on SLURM job configuration."
  elif [[ $DEFAULT_THREAD_COUNT ]]
  then
    argThreads=$DEFAULT_THREAD_COUNT
    echo "Setting --threads=$argThreads based on DEFAULT_THREAD_COUNT"
  else
    argThreads=1
  fi
fi

# Start building a string with all the additional cufflinks params
# (The only things not included here are the input/output paths)
cuffParams="-p $argThreads"

#  --gff - If not specified, check for FLY_KNOWN_GFF, or throw an error
#		Specified file must exist
if [[ ! $argGFF ]]
then
  if [[ $FLY_KNOWN_GFF ]]
  then
    argGFF=$FLY_KNOWN_GFF
    echo "Inferred --gff=$argGFF from FLY_KNOWN_GFF"
  elif [[ $argMode != "denovo" ]]
  then
    echo "ERROR: Must specify GFF file with --gff parameter or FLY_KNOWN_GFF environment variable"
    exit 1
  fi
fi

if [[ $argGFF ]]
then
  if [[ ! -e $argGFF ]]
  then
    echo "ERROR: Missing required file $argGFF"
    exit 1
  fi
fi

# --mode = Must be one of: guide, known, denovo
#		If not specified, default to guide (use -g $gff param)
if [[ ! $argMode ]]
then
  argMode="guide"
  echo "Setting --mode=$argMode by default."
fi

# Setup appropriate cufflinks param based on --mode and --gff
if [[ $argMode == "guide" ]]
then
  echo "Using $argGFF as a guide transcriptome, but allowing for novel transcript discovery (RABT algorithm)."
  cuffParams="$cuffParams -g $argGFF"
elif [[ $argMode == "known" ]]
then
  echo "Quantitating transcriptome against reference assembly in $argGFF, no novel transcripts allowed."
  cuffParams="$cuffParams -G $argGFF"
elif [[ $argMode == "denovo" ]]
then
  echo "Performing de novo transcript assembly - no guide transcriptome."
  if [[ $argGFF ]]
  then
  	echo "WARNING: --gff=$argGFF is suppressed due to --mode=$argMode option."
  fi
else
  echo "ERROR: --mode=$argMode is not a valid option."
  exit 1
fi

# --library-type - default to fr-firststrand, must be one of [fr|ff]-[unstranded|firststrand|secondstrand]
# TO DO: Create a project-wide variable for this?
if [[ ! $argLibraryType ]]
then
  argLibraryType="fr-firststrand"
  echo "Setting --library-type=$argLibraryType by default"
fi

ltPrefix=${argLibraryType:0:3}
ltSuffix=`echo "$argLibraryType" | sed 's/^f[fr]-//g'`

if [[ ("$ltPrefix" != "fr-") && ("$ltPrefix" != "ff-") ]]
then
  echo "ERROR: --library-type=$argLibraryType is not valid, see cufflinks documentation"
  exit 1
fi

if [[ ("$ltSuffix" != "unstranded") && ("$ltSuffix" != "firststrand") && ("$ltSuffix" != "secondstrand") ]]
then
  echo "ERROR: --library-type=$argLibraryType is not valid, see cufflinks documentation"
  exit 1
fi

cuffParams="$cuffParams --library-type $argLibraryType"

# --max-bundle-frags - default is 1000000, Just make sure it's an integer > 0, warn if less than default
# TO DO: Create a project-wide variable for this?
# Might also be wise to automatically scale this to total reads...

# NOTE: If not specified, this now gets left off the cuffParams string - so goes with cufflinks default
if [[ $argMaxBundleFrags ]]
then
  if [[ ! ("$argMaxBundleFrags" =~ ^[0-9]+$ ) ]]
  then
    echo "ERROR: --max-bundle-frags must be a positive integer"
    exit 1
  fi
  
  if [[ "$argMaxBundleFrags" -lt 1000000 ]]
  then
    echo ""
    echo "WARNING: --max-bundle-frags=$argMaxBundleFrags is lower than default value, may exclude high expression genes!"
    echo ""
  fi
  
  cuffParams="$cuffParams --max-bundle-frags $argMaxBundleFrags"
fi

# --max-bundle-length - default is 3500000, Just make sure it's an integer > 0, warn if less than default
# TO DO: Create a project-wide variable for this?
# Might also be wise to automatically scale this to total reads...

# NOTE: If not specified, this now gets left off the cuffParams string - so defers to cufflinks default
if [[ $argMaxBundleLength ]]
then
  if [[ ! ("$argMaxBundleLength" =~ ^[0-9]+$ ) ]]
  then
    echo "ERROR: --max-bundle-length must be a positive integer"
    exit 1
  fi
  
  if [[ "$argMaxBundleLength" -lt 3500000 ]]
  then
    echo ""
    echo "WARNING: --max-bundle-length=$argMaxBundleLength is lower than default value, may exclude gene-dense regions!"
    echo ""
  fi
  
  cuffParams="$cuffParams --max-bundle-length $argMaxBundleLength"
fi

# --extraparams - if specified, just append to the cuffParams string, no error checking
if [[ $argExtraParams ]]
then
  cuffParams="$cuffParams $argExtraParams"
fi


# --dryrun - Just report if this is the case
if [[ $argDryRun ]]
then
  echo "---- DRY RUN ----"
  echo "Tasks will not actually run, there will be no output!"
else
  mkdir -p $argOutput
fi


# TO DO: Check for dependencies and output the versions, e.g.
dependencies=("cufflinks")
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
    elif [[ $tool == "cufflinks" ]]
    then
      $tool 2>&1 | head -n 2
    else
      $tool --version
    fi
    echo ""
  fi
done


# -- MAIN TASK -- #

# Build and run the main cufflinks command
# Prints the cmd to run, then runs it (output from cmd goes to stdout/stderr)

# Generate full path to aligned BAM file
alignedBam=$argInput/$argSample/$argBamFile

# Make sure the file exists
if [[ ! -e $alignedBam ]]
then
  echo "ERROR: Missing $alignedBam"
  exit 1
else
  echo "Assembling aligned reads from $alignedBam"
fi

# Specify Output path, check if results exist already
outputPath=$argOutput/$argSample
if [[ -e $outputPath/transcripts.gtf ]]
then
  echo "ERROR: Cufflinks output file already exists:"
  echo "$outputPath/transcripts.gtf"
  echo "Remove or rename this file before running this task"
  exit 1
else
  echo "Cufflinks output will be in: $outputPath/"
fi
  
# Because the output log is very long, make sure to pipe command output to separate file
logFile=$outputPath/cufflinks.log

echo ""
echo "RUNNING: cufflinks --output-dir $outputPath $cuffParams $alignedBam &> $logFile"
echo "Check $logFile for program status"

if [[ ! $argDryRun ]]
then
  echo -ne "Started at:\t"
  date
  echo ""
  
  mkdir -p $outputPath
  
  cufflinks --output-dir $outputPath $cuffParams $alignedBam &> $logFile
  cufflinks_exit_status=$?
  
  echo -ne "Finished at:\t"
  date
  echo ""
  
  warnFile=$outputPath/cufflinks_warnings.log
  grep 'Warning' $logFile > $warnFile
  
  echo "Check $warnFile for warning messages."
  
  if [[ "cufflinks_exit_status" -ne "0" ]]
  then
    echo "ERROR: cufflinks returned non-zero exit status = $cufflinks_exit_status"
    exit $cufflinks_exit_status
  fi
fi

echo "Entire script completed successfully."
