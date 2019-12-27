#!/bin/bash

# LJE - 3/24/15
#
# run_batch_task.sh - Script for running any generic task over a batch of a samples, with automated sbatch submission
#
# To run on a single sample:
#  ./run_batch_task.sh --batch=batchname [--sample=X|--line=X] --task='CMDSTRING' [options...]
# To run on all samples in serial (no cluster submission):
#  ./run_batch_task.sh --batch=batchname --task='CMDSTRING' [options...]
# To submit batch of jobs to cluster:
#  sbatch --array=1-L [-c N] run_batch_task.sh --batch=batchname --task='CMDSTRING' [options...]
# OR:
#  ./run_batch_task.sh --batch=batchname --task='CMDSTRING' --array [--sbatchopt='-c N'] [options...]
# 
# batchname.txt must be a tab-delimited file with one line per sequence file
# with columns: SampleID, SampleName, FlowCellID, Lane, Index, Barcode, RunDate
# In general, SampleID is the unique identifier used to name each of the processed files and matches up to the original fastq file
#
# CMDSTRING must specify the task to run. Prior to executing this command, 
# the following variables will be defined based on the current line of batchname.txt:
# $sampleID, $sampleName, $flowcell, $lane, $index, $barcode, $rundate
# These variables can be used in specifying the CMDSTRING, just make sure to protect it with single quotes
# NOTE: batchname is also be used for output paths, it can be passed to the primary task as '$argBatch' variable
#
# Generic Options:
#  -b --batch=batchname - This should refer to a tab-delimited text file with sample info (first column must be unique ID)
#		The file must end in .txt and the extension can be excluded when specifying this parameter
#		or will be dropped automatically when assigning names to batch subdirectories
#  -t --task='CMDSTRING' - This will be run through eval after identifying the appropriate set of files corresponding to sample ID
#		Within CMDSTRING, use "${mySampleFiles[@]}" to insert the list of sample files
#		and use "mySampleName" to insert the sample name
#		Always use '' to protect the CMDSTRING here so it doesn't get evaluated until later
#  -p --paired - Indicates that this batch is paired-end data, this creates a second variable, $pairedID, for use in the task command
#  -a --array - Indicates that this instance should initiate the slurm array command itself
#  --sbatchopt - String of parameters to pass to sbatch command, best to protect with '' if specified
#  -s --sample=sampleID - Must match the first field in a line from batchname.txt, this will run the task on just that sample
#  -l --line=X - Alternate to sampleID, runs the task on just the sample in line X of batchname.txt
#  -c --threads=N - Specifies the number of threads to use, should also be passed to -c param of sbatch
#		If unspecified here, script checks $SLURM_JOB_CPUS_PER_NODE, then $DEFAULT_THREAD_COUNT, or defaults to 1
#		The value here can be passed to primary task as $argThreads, however it's better to make the primary task script
#		capable of auto-sensing $SLURM_JOB_CPUS_PER_NODE
#  -d --dryrun - Parameter indicating that task should not actually be run, but any pre-task checks should be run,
#		Note: this will not run primary task at all, just show how the final execution string will look.
#

# parse command line options
args=`getopt -o "b:t:s:l:c:apd" -l "batch:,task:,paired,array,sbatchopt:,sample:,line:,threads:,dryrun" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -b|--batch)
      argBatch=$2
      shift 2;;
    
    -t|--task)
      argTask=$2
      shift 2;;
    
    -p|--paired)
      argPaired=1
      shift 1;;
    
    -a|--array)
      argArray=1
      shift 1;;
    
    --sbatchopt)
      argSbatchOpt=$2
      shift 2;;
    
    -s|--sample)
      argSample=$2
      shift 2;;
    
    -l|--line)
      argLine=$2
      shift 2;;
    
    -c|--threads)
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

# Attempt to load project info - contains defaults for many variables
if [[ -e project_info.sh ]]
then
  source project_info.sh
  echo "Loaded project_info.sh"
else
  echo "No project_info.sh, attempting to proceed without it..."
fi

# Check each parameter
# Specify required params here and set defaults for others

# --batch - Required, must point to an existing file, if full file name is given, drop the ".txt" extension
if [[ ! $argBatch ]]
then
  echo "ERROR: Must specify a batch name!"
  exit 1
fi

# Remove trailing .txt if specified
argBatch=`echo "$argBatch" | sed 's/\.txt$//g'`

if [[ ! -e $argBatch.txt ]]
then
  echo "ERROR: Cannot find $argBatch.txt"
  exit 1
fi

# Needed for both range checking and automated array submission
batchSize=`wc -l $argBatch.txt | awk '{print $1}'`

# --task - REQUIRED
if [[ ! $argTask ]]
then
  echo "ERROR: Must define a --task parameter"
  exit 1
fi

# --sample - Must match first field on exactly 1 line in $argBatch.txt
if [[ $argSample ]]
then
  sampleMatch=`awk -F'\t' '($1=="'$argSample'")' $argBatch.txt | wc -l`
  if [[ "$sampleMatch" -lt "1" ]]
  then
    echo "ERROR: $argSample does not match any entries in $argBatch.txt"
    exit 1
  elif [[ "$sampleMatch" -gt "1" ]]
  then
    echo "ERROR: $argSample matched multiple entries in $argBatch.txt"
    exit 1
  fi
fi

# --line - Infer from $SLURM_ARRAY_TASK_ID if specified, Must be between 1 and the number of lines in $argBatch.txt
# If --line is specified or inferred from SLURM, there should not be a --sample parameter

if [[ $SLURM_ARRAY_TASK_ID ]]
then
  if [[ $argSample || $argLine ]]
  then
    echo "ERROR: Do not use the --sample or --line parameters when running in array mode"
    exit 1
  else
    argLine=$SLURM_ARRAY_TASK_ID
    echo "Inferring --line=$argLine from SLURM_ARRAY_TASK_ID"
  fi
fi

if [[ $argLine ]]
then
  if [[ $argSample ]]
  then
    echo "ERROR: Do not specify both --sample and --line parameters"
    exit 1
  fi
  if [[ "$argLine" -lt "1" ]]
  then
    echo "ERROR: --line must specify a positive number"
    exit 1
  elif [[ "$argLine" -gt "$batchSize" ]]
  then
    echo "ERROR: cannot process --line=$argLine, $argBatch.txt only has $batchSize entries"
    exit 1
  fi
fi

# --array - Not compatible with --line, --name, or $SLURM_JOB_ID
if [[ $argArray ]]
then
  if [[ $argIndex || $argName || $SLURM_JOB_ID ]]
  then
    echo "ERROR: Cannot use --array parameter with --index, --name, or when passing job directly to SLURM controller"
    echo "The purpose of the --array option is to automatically calculate the array range and call sbatch within this shell"
    exit 1
  fi
fi 

# --threads - If not specified, check for SLURM_JOB_CPUS_PER_NODE or DEFAULT_THREAD_COUNT, or default to 1
if [[ ! $argThreads ]]
then
  if [[ $SLURM_JOB_CPUS_PER_NODE ]]
  then
    argThreads=$SLURM_JOB_CPUS_PER_NODE
  elif [[ $DEFAULT_THREAD_COUNT ]]
  then
    argThreads=$DEFAULT_THREAD_COUNT
  else
    argThreads=1
  fi
fi

# --dryrun - Just report if this is the case
if [[ $argDryRun ]]
then
  echo "---- DRY RUN ----"
  echo "Tasks will not actually run, there will be no output!"
  ech "----- DRY RUN ----"
  echo ""
fi


# -- CORE SUBROUTINES -- #

# Define the core function to parse lines from batch file and build final task command
# Takes one argument: a single line from batch file (NOTE: Must be quoted)
# The fields are parsed and stored in variables accessible to the command string in $argTask
# Prints the final command to run, then runs it (output goes to stdout/stderr)
run_task () {
  line=$1
  sampleID=`echo "$line" | cut -f 1`
  sampleName=`echo "$line" | cut -f 2`
  flowcell=`echo "$line" | cut -f 3`
  lane=`echo "$line" | cut -f 4`
  index=`echo "$line" | cut -f 5`
  barcode=`echo "$line" | cut -f 6`
  rundate=`echo "$line" | cut -f 7`
  
  # If paired data, replace "_R1_" with "_R2_" to get pairedID
  # NOTE: Reversing before sed guarantees that the LAST instance of _R1_ gets replace with _R2_ (in case _R1_ also part of sample name)
  pairedID=`echo $sampleID | rev | sed 's/[_]1R[_]/_2R_/' | rev`
  
  # Check that sampleID matches other fields
  # This is really just to validate my change in batch file format on 11/17/14
  # But probably a good idea to keep this formatting rule going forward
  # NO LONGER ENFORCING THIS FORMAT
  # if [[ $sampleID != "${sampleName}_${barcode}_L00${lane}_${flowcell}" ]]
  # then
  #   echo "ERROR: sampleID \"$sampleID\" does not match other fields"
  #   echo "Expected \"${sampleName}_${barcode}_L00${lane}_${flowcell}\""
  #   exit 1
  # fi
  
  eval echo \"RUNNING: $argTask\"
  echo ""
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
    
    eval $argTask
    task_exit_code=$?
        
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
    
    if [[ "$task_exit_code" -ne "0" ]]
    then
      echo "ERROR: Primary command returned non-zero exit status = $task_exit_code"
      exit $task_exit_code
    fi
  fi
}


# -- MAIN -- #

# If --array param was set, then just rebuild an sbatch command around the params given here
if [[ $argArray ]]
then
  # Array range: 1-${#sampleNames[@]}
  sbatchCmd="sbatch --array=1-$batchSize"
  if [[ $argSbatchOpt ]]
  then
    sbatchCmd+=" $argSbatchOpt"
  fi
  # The command itself is to call this script with same params except for --array and --sbatchopt
  sbatchCmd+=" run_batch_task.sh"
  if [[ $argBatch ]]
  then
    sbatchCmd+=" --batch=$argBatch"
  fi
  if [[ $argTask ]]
  then
    sbatchCmd+=" --task='$argTask'"
  fi
  if [[ $argSample ]]
  then
    sbatchCmd+=" --sample=$argSample"
  fi
  if [[ $argLine ]]
  then
    sbatchCmd+=" --line=$argLine"
  fi
  if [[ $argThreads ]]
  then    
    sbatchCmd+=" --threads=$argThreads"
  fi
  if [[ $argDryRun ]]
  then
    sbatchCmd+=" --dryrun"
  fi
  
  echo "Submitting job array to SLURM:"
  echo "$sbatchCmd"
  
  eval $sbatchCmd
  sbatch_exit_code=$?
  
  exit $sbatch_exit_code
fi

# Otherwise proceed to running the actual task  

# Report where/when we are
CURWD=`pwd`
CURHOST=`hostname`
echo "Running in directory $CURWD on $CURHOST"

# Check job manager variables
if [[ $SLURM_JOBID ]]
then
  echo "Run through SLURM, Job ID = $SLURM_JOBID"
fi
if [[ $SLURM_ARRAY_JOB_ID ]]
then
  echo "Run as part of a SLURM job array: Array Job ID = $SLURM_ARRAY_JOB_ID, Task ID = $SLURM_ARRAY_TASK_ID"
fi

# Figure out what mode to run in:
if [[ $argSample || $argLine ]]
then
  # Run in sample-specific mode
  if [[ $argSample ]]
  then
    batchline=`awk -F'\t' '($1 == "'$argSample'")' $argBatch.txt`
  elif [[ $argLine ]]
  then
    batchline=`head -n $argLine $argBatch.txt | tail -n 1`
  fi
  
  # Call core task function
  run_task "$batchline"
else
  # Run in loop mode - process all lines in $argBatch.txt
  echo "Processing all samples in $argBatch.txt in serial mode"
  while read batchline
  do
    run_task "$batchline"
  done < $argBatch.txt
fi

echo "run_batch_task.sh completed successfully."
