#!/bin/bash
# NOTE: Do not set -e here, or will crash when checking samtools version, haven't found a good work around for this yet

# LJE - 4/2/15
#
# run_cuffmerge.sh - Run cuffmerge on a batch of cufflinks results
# By its nature, this only requires a single job for a whole batch, do NOT run with run_batch_task.sh, 
# just queue up directly with sbatch if running on the cluster
#
# Usage:
# ./run_cuffmerge.sh --batch=batchname [options]
#
# Options:
#  --input=path - Path to the input directory, default is $ASSEMBLED_PATH, if defined in project_info.pl, otherwise defaults to cufflinks
#  --filename=string - The standard gtf file name in the cufflinks directory, defaults to "transcripts.gtf",
#		Change this only if some pre-filtering has been done on theses cufflinks output files prior to merging
#  --output=path - Path to the output directory, default is $CUFFMERGE_PATH, as defined in project_info.sh (if missing, defaults to "cuffmerge"
#		The actual output goes in a subdirectory named sampleID in the path indicated here
#		If STAR output for this sample exists already, this script will quit with an error
#  --batch=batchname - Convenience option, adds a matching subdirectory to both the input and output paths
#  --gff=path/to/gff - File to use for gene annotations, default is taken from FLY_KNOWN_GFF, 
#			NOTE: Currently only using this for cuffcompare run, NOT including it in cuffmerge step
#			TO DO: Create a parameter to skip the cuffcompare step?
#  --prevgtf - (Optional) Specify an extra GTF file of previous transcript models (NOT part of the "Known" GFF assembly)
#			This was specifically used to include Wen's previous NTRs in the DGRP Baseline assembly,
#			But it could be used to include previous novel NTRs in NTR assemblies for additional conditions
#  --genome=path/to/fasta - File to use for genome sequence, default is taken from FLY_GENOME_FASTA
#			Causes an error if that value is missing
#  --extraparams - Pass additional parameters to cufflinks that are aren't currently part of the options here
#			Just make sure to protect with single quotes
#  --threads=N - Specifies the number of threads to use, if unspecified, will check SLURM environment, then DEFAULT_THREAD_COUNT, otherwise defaults to 1
#  --dryrun - Parameter indicating that task should not actually be run, but any pre-task checks should be run,
#		and the final formatted command should be echo'd. This is useful for testing and debugging.
#		Basically, if dryrun is specified, there should be no output or file manipulation (code below must honor this though)
#

# -- COMMAND LINE PARAMETERS -- #

# parse command line options
# NOTE: all variables specified here start with "arg"
args=`getopt -o "i:f:o:b:m:t:d" -l "input:,filename:,output:,batch:,gff:,prevgtf:,mode:,library-type:,max-bundle-frags:,max-bundle-length:,extraparams:,threads:,dryrun" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -i|--input)
      argInput=$2
      shift 2;;
    
    -f|--filename)
      argFilename=$2
      shift 2;;
    
    -o|--output)
      argOutput=$2
      shift 2;;
        
    -b|--batch)
      argBatch=$2
      shift 2;;
    
    --gff)
      argGFF=$2
      shift 2;;
    
    --prevgtf)
      argPrevGTF=$2
      shift 2;;
    
    --genome)
      argGenome=$2
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

#  --input=path - Path to the input directory, default is $ASSEMBLED_PATH, if defined in project_info.pl, otherwise defaults to cufflinks
#			Directory must exist
if [[ ! $argInput ]]
then
  if [[ $ASSEMBLED_PATH ]] 
  then
    argInput=$ASSEMBLED_PATH
    echo "Inferring --input=$argInput based on ASSEMBLED_PATH in project_info.sh"
  else
    argInput="cufflinks"
    echo "Setting --input=$argInput by default"
  fi
fi

#  --filename=string - The standard gtf file name in the cufflinks directory, defaults to "transcripts.gtf",
#		Change this only if some pre-filtering has been done on theses cufflinks output files prior to merging
if [[ ! $argFilename ]]
then
  argFilename="transcripts.gtf"
  echo "Setting --filename=$argFilename by default"
fi

# --output, if not specified check for CUFFMERGE_PATH, or default to "cuffmerge"
if [[ ! $argOutput ]]
then
  if [[ $CUFFMERGE_PATH ]]
  then
    argOutput=$CUFFMERGE_PATH
    echo "Inferring --output=$argOutput based on CUFFMERGE_PATH in project_info.sh"
  else
    argOutput=cuffmerge
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

#  --gff - If not specified, check for FLY_KNOWN_GFF, or throw an error
#		Specified file must exist
# NOTE: This is only used for cuffcompare, NOT cuffmerge!
if [[ ! $argGFF ]]
then
  if [[ $FLY_KNOWN_GFF ]]
  then
    argGFF=$FLY_KNOWN_GFF
    echo "Inferred --gff=$argGFF from FLY_KNOWN_GFF"
  else
    echo "ERROR: Must specify GFF file with --gff parameter or FLY_KNOWN_GFF in project_info.pl"
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

#  --genome=path/to/fasta - File to use for genome sequence, default is taken from FLY_GENOME_FASTA
#			Causes an error if that value is missing
if [[ ! $argGenome ]]
then
  if [[ $FLY_GENOME_FASTA ]]
  then
    argGenome=$FLY_GENOME_FASTA
    echo "Inferred --genome=$argGenome from FLY_GENOME_FASTA"
  else
    echo "ERROR: Must specify genome fasta file with --genome param or FLY_GENOME_FASTA in project_info.pl"
    exit 1
  fi
fi

if [[ ! -e $argGenome ]]
then
  echo "ERROR: Missing required file $argGenome"
  exit 1
fi

# Start the params string:
cuffParams="-o $argOutput -s $argGenome -p $argThreads --keep-tmp"

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
dependencies=("cuffmerge" "cuffcompare")
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
    elif [[ $tool == "cufflinks" || $tool == "cuffcompare" ]]
    then
      $tool 2>&1 | head -n 2
    else
      $tool --version
    fi
    echo ""
  fi
done


# -- MAIN TASK -- #

# Create the output directory
if [[ ! $argDryRun ]]
then
  mkdir -p $argOutput
fi



# Make sure none of the output files exist yet
gtfList=$argOutput/cuffmerge_gtf_list.txt
mergeFile=$argOutput/merged.gtf
compFile=$argOutput/cuffcompare.combined.gtf

if [[ -e $gtfList ]]
then
  echo "ERROR: Cuffmerge input list already exists:"
  echo "$gtfList"
  echo "Remove or rename this file before running this task"
  exit 1
fi

if [[ -e $mergeFile ]]
then
  echo "ERROR: Cuffmerge output file already exists:"
  echo "$mergeFile"
  echo "Remove or rename this file before running this task"
  exit 1
fi

if [[ -e $compFile ]]
then
  echo "ERROR: Cuffcompare output file already exists:"
  echo "$compFile"
  echo "Remove or rename this file before running this task"
  exit 1
fi

# Generate the cuffmerge_gtf_list.txt file
echo "Generating list of GTF files to merge:"
echo "RUNNING: ls $argInput/*/$argFilename > $gtfList"

if [[ ! $argDryRun ]]
then
  ls $argInput/*/$argFilename > $gtfList
  
  if [[ ! -e $gtfList ]]
  then
    echo "ERROR: Was not able to generate $gtfList"
    exit 1
  fi
  
  # Spike in previous GTF file if specified
  if [[ $argPrevGTF ]]
  then
    echo "Adding $argPrevGTF to end of $gtfList"
    echo "$argPrevGTF" >> $gtfList
  fi
  
fi

# Run cuffmerge step

# Build and run the main cuffmerge command
# Prints the cmd to run, then runs it (output from cmd goes to stdout/stderr)

# Because the output log is very long, make sure to pipe command output to separate file
logFile=$argOutput/cuffmerge.log

echo ""
echo "RUNNING: cuffmerge $cuffParams $gtfList &> $logFile"
echo "Check $logFile for program status"

if [[ ! $argDryRun ]]
then
  echo -ne "Started at:\t"
  date
  echo ""
  
  mkdir -p $argOutput
  
  cuffmerge $cuffParams $gtfList &> $logFile
  cuffmerge_exit_status=$?
  
  echo -ne "Finished at:\t"
  date
  echo ""
  
  warnFile=$argOutput/cuffmerge_warnings.log
  grep -f <(echo -e "Warning\nError") -i $logFile > $warnFile
  
  echo "Check $warnFile for warning messages."
  
  if [[ "cuffmerge_exit_status" -ne "0" ]]
  then
    echo "ERROR: cuffmerge returned non-zero exit status = $cuffmerge_exit_status"
    exit $cuffmerge_exit_status
  fi
  
  if [[ ! -e $mergeFile ]]
  then
    echo "ERROR: Missing output from cuffmerge:"
    echo "$mergeFile"
    echo "Check $logFile for more details"
    exit 1
  fi
  
  # Clean up the $argOutput/tmp/ directory
  if [[ -e $argOuput/tmp/ ]]
  then
    echo "Cleaning up: $argOutput/tmp/"
    rm -R $argOutput/tmp/
  fi
fi

# Run cuffcompare step
echo ""
echo "RUNNING: cuffcompare -r $argGFF -s $argGenome -o $argOutput/cuffcompare $mergeFile"

if [[ ! $argDryRun ]]
then
  echo -ne "Started at:\t"
  date
  echo ""
  
  cuffcompare -r $argGFF -s $argGenome -o $argOutput/cuffcompare $mergeFile
  cuffcompare_exit_status=$?
  
  echo -ne "Finished at:\t"
  date
  echo ""
  
  if [[ "$cuffcompare_exit_status" -ne "0" ]]
  then
    echo "ERROR: cuffcompare returned non-zero exit status = $cuffcompare_exit_status"
    exit $cuffcompare_exit_status
  fi
  
  if [[ ! -e $compFile ]]
  then
    echo "ERROR: Missing output from cuffcompare:"
    echo "$compFile"
    exit 1
  fi
fi

echo "Entire script completed successfully."
