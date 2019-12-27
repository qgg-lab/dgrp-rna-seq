#!/bin/bash

# LJE - 3/24/15
# Major Update - 12/11/15 - Now has modes for different filtering tasks
#
# run_bwa.sh - Script for running BWA on sequence files for DGRP Baseline RNA-seq project
# This was originally intended for rRNA filtering, so it uses the rRNA genome by default
# And automatically separates the mapped and unmapped reads back into fastq files
# This version does NOT read batch files or handle SLURM job arrays - use with run_batch_task.sh for that functionality
#
# Main update: this script is now also to be used for microbe filtering (after ribosomal filtering)
# This can already be done by modifying multiple parameters, 
# so I am simply adding a "--mode" parameter that alters the appropriate defaults.
# This can be extended in the future to other preset filtering tasks.
#
# Usage:
# ./run_bwa.sh [options] --sample=sampleID
# 
# Options:
#  --sample=sampleID - Refers to the file to process
#  --batch=batchname - If specified, gets appended to end of both --input and --output paths
#  --mode=[ribo|microbe|repeat] - Adjusts default values for other parameters, default mode is ribo (ribosomal filtering)
#		- Set to "microbe" for microbial filtering step
#		- Set to "repeat" for repeat filtering step
#  --input=path - Path to the input directory, defaults depend on mode:
#		ribo: defaults to $TRIMMED_FASTQ_PATH, $PROC_FASTQ_PATH, or trimmed/ if none are defined
#		microbe: defaults to $RIBO_FILTERED_PATH, or ribofiltered/ if neither are defined
#		repeat: defaults $MICROBE_FILTERED_PATH, or microbefiltered/ if neither are defined
# 		This directory must contain a file called (sampleID)_(insuffix).fastq(.gz) (or paired sequence files if --paired)
#  --insuffix=string - String to use as suffix of of input file (after sampleID, do not include .gz for compressed files)
#		Default (ribo mode): _trimmed.fastq
#		Default (microbe mode): _filtered.fastq
#		Default (repeat mode): _filtered.fastq
#  --output=path - Path to the output directory, default depends on --mode:
#		ribo mode: Defaults to $BWA_RIBO_PATH, $BWA_PATH, or bwa_rRNA/ if neither are defined
#       microbe mode: Defaults to $BWA_MICROBE_PATH, or bwa_microbe/
#		repeat mode: Defaults to $BWA_REPEAT_PATH, or bwa_repeat/
#		Final output goes in outpath/sampleID/ - if BWA output exists here already, this script will not run
#  --genome=path - Path to BWA genome index directory. Default depends on --mode:
#		ribo: checks $FLY_RIBOSOMAL_RNA, then hits an error...
#		microbe: checks $MICROBE_DB, then hits an error...
#		repeat: checks $REPEAT_DB, then hits an error
#  --paired - Indicates paired read data as input
#  --mapped=path - Path to store the mapped reads back into fastq format as (sampleID)(_1|2)_mapped.fastq.gz, defaults to the main output directory
#  --unmapped=path - Path to store the unmapped reads back into fastq format as (sampleID)(_1|2)_filtered.fastq.gz, 
#		ribo mode: defaults to $RIBO_FILTERED_PATH if specified in project_info.sh
#		microbe mode: defaults to $MICROBE_FILTERED_PATH if specified in project_info.sh
#		repeat mode: defaults to $REPEAT_FILTERED_PATH if specified in project_info.sh
#		In either case, if no default in project_info.sh, defaults to the main output directory
#  --threads=N - Specifies the number of threads to use, should also be passed to -c param of sbatch
#		If unspecified here, script checks DEFAULT_THREAD_COUNT, or defaults to 1
#		TO DO: Check SLURM environment variable instead of DEFAULT_THREAD_COUNT
#  --dryrun - Parameter indicating that task should not actually be run, but any pre-task checks should be run,
#		and the final formatted command should be echo'd. This is useful for testing and debugging.
#		Basically, if dryrun is specified, there should be no output or file manipulation (code below must honor this though)
#

# parse command line options
args=`getopt -o "b:s:l:t:i:o:g:pm:u:d" -l "batch:,sample:,mode:,line:,threads:,input:,insuffix:,output:,genome:,paired,mapped:,unmapped:,dryrun" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -s|--sample)
      argSample=$2
      shift 2;;
    
    --mode)
      argMode=$2
      shift 2;;
    
    -i|--input)
      argInput=$2
      shift 2;;
    
    --insuffix)
      argInSuffix=$2
      shift 2;;
    
    -o|--output)
      argOutput=$2
      shift 2;;
    
    -b|--batch)
      argBatch=$2
      shift 2;;
    
    -g|--genome)
      argGenome=$2
      shift 2;;
    
    -p|--paired)
      argPaired=1
      shift 1;;
    
    -m|--mapped)
      argMapped=$2
      shift 2;;
    
    -u|--unmapped)
      argUnmapped=$2
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

# --sample - Required
if [[ ! $argSample ]]
then
  echo "ERROR: Must specify --sample parameter!"
  exit 1
fi

# --mode - Default to ribo if not set, otherwise make sure it's a valid value
if [[ ! $argMode ]]
then
  argMode="ribo"
  echo "Setting --mode=$argMode (rRNA filtering step) by default"
elif [[ ($argMode == "ribo") || ($argMode == "microbe") || ($argMode == "repeat") ]]
then
  echo "Setting parameter defaults for $argMode mode"
else
  echo "ERROR: --mode=$argMode is not a valid mode, options are: ribo, microbe, repeat"
  exit 1
fi

# --input, if not specified check for appropriate default depending on mode
if [[ ! $argInput ]]
then
  if [[ $argMode == "microbe" ]]
  then
    # microbe mode default: $RIBO_FILTERED_PATH, or ribofiltered/
    if [[ $RIBO_FILTERED_PATH ]]
    then
      argInput=$RIBO_FILTERED_PATH
      echo "Inferring --input=$argInput based on RIBO_FILTERED_PATH in project_info.sh"
    else
      argInput="ribofiltered"
      echo "Setting --input=$argInput by default"
    fi
  elif [[ $argMode == "repeat" ]]
  then
    # repeat mode default: $MICROBE_FILTERED_PATH, or microbefiltered/
    if [[ $MICROBE_FILTERED_PATH ]]
    then
      argInput=$MICROBE_FILTERED_PATH
      echo "Inferring --input=$argInput based on MICROBE_FILTERED_PATH in project_info.sh"
    else
      argInput="microbefiltered"
      echo "Setting --input=$argInput by default"
    fi
  else
  	# ribo mode settings by default
  	if [[ $TRIMMED_FASTQ_PATH ]]
  	then
  	  argInput=$TRIMMED_FASTQ_PATH
  	  echo "Inferring --input=$argInput based on TRIMMED_FASTQ_PATH in project_info.sh"
  	elif [[ $PROC_FASTQ_PATH ]] 
  	then
      argInput=$PROC_FASTQ_PATH
      echo "Inferring --input=$argInput based on PROC_FASTQ_PATH in project_info.sh"
  	else
      argInput=trimmed
      echo "Setting --input=$argInput by default"
  	fi
  fi
fi

# --insuffix, if not specified, default to _trimmed.fastq (ribo mode) or _filtered.fastq (other modes)
if [[ ! $argInSuffix ]]
then
  if [[ ($argMode == "microbe") || ($argMode == "repeat") ]]
  then
    argInSuffix="_filtered.fastq"
    echo "Setting --insuffix=$argInSuffix by default ($argMode mode)"
  else
    argInSuffix="_trimmed.fastq"
    echo "Setting --insuffix=$argInSuffix by default (ribo mode)"
  fi
fi

# --output, if not specified set defaults based on --mode
#		ribo mode: Defaults to $BWA_RIBO_PATH, $BWA_PATH, or bwa_rRNA/ if neither are defined
#       microbe mode: Defaults to $BWA_MICROBE_PATH, or bwa_microbe/
#		repeat mode: Defaults to $BWA_REPEAT_PATH, or bwa_repeat/
#		Final output goes in outpath/sampleID/ - if BWA output exists here already, this script will not run
#
if [[ ! $argOutput ]]
then
  if [[ $argMode == "microbe" ]]
  then
    if [[ $BWA_MICROBE_PATH ]]
    then
      argOutput=$BWA_MICROBE_PATH
      echo "Inferring --output=$argOutput based on BWA_MICROBE_PATH in project_info.sh"
    else
      argOutput="bwa_microbe"
      echo "Setting --output=$argOutput by default (microbe mode)"
    fi
  elif [[ $argMode == "repeat" ]]
  then
    if [[ $BWA_REPEAT_PATH ]]
    then
      argOutput=$BWA_REPEAT_PATH
      echo "Inferring --output=$argOutput based on BWA_REPEAT_PATH in project_info.sh"
    else
      argOutput="bwa_repeat"
      echo "Setting --output=$argOutput by default (repeat mode)"
    fi
  else
    if [[ $BWA_RIBO_PATH ]]
    then
      argOutput=$BWA_RIBO_PATH
      echo "Inferring --output=$argOutput based on BWA_RIBO_PATH in project_info.sh"
    elif [[ $BWA_PATH ]]
    then
      argOutput=$BWA_PATH
      echo "Inferring --output=$argOutput based on BWA_PATH in project_info.sh"
    else
      argOutput="bwa_rRNA"
      echo "Setting --output=$argInput by default (ribo mode)"
    fi
  fi
fi

# --batch - Optional - append to the end of argInput and argOutput...
if [[ $argBatch ]]
then
  argInput=$argInput/$argBatch
  argOutput=$argOutput/$argBatch
  echo "Using matching batch subdirectories for input and output, final paths:"
  echo "Input: $argInput"
  echo "Output: $argOutput"
fi

# Make sure argInput is valid
if [[ ! -d $argInput ]]
then
  echo "ERROR: Input directory does not exist"
  echo "MISSING: $argInput"
  exit 1
fi

# --genome=path, if not specified set defaults based on mode
if [[ ! $argGenome ]]
then
  if [[ $argMode == "microbe" ]]
  then
    if [[ $MICROBE_DB ]]
    then
      argGenome=$MICROBE_DB
      echo "Inferring --genome=$argGenome from MICROBE_DB in project_info.sh"
    else
      echo "ERROR: Must specify genome index with --genome, or in project_info.sh with MICROBE_DB (for --mode=microbe)"
      exit 1
    fi
  elif [[ $argMode == "repeat" ]]
  then
    if [[ $REPEAT_DB ]]
    then
      argGenome=$REPEAT_DB
      echo "Inferring --genome=$argGenome from REPEAT_DB in project_info.sh"
    else
      echo "ERROR: Must specify genome index with --genome, or in project_info.sh with REPEAT_DB (for --mode=repeat)"
      exit 1
    fi
  else
    # ribo mode defaults
    if [[ $FLY_RIBOSOMAL_RNA ]]
    then
      argGenome=$FLY_RIBOSOMAL_RNA
      echo "Inferring --genome=$argGenome from FLY_RIBOSOMAL_RNA in project_info.sh"
    else
      echo "ERROR: Must specify genome index with --genome, or in project_info.sh with FLY_RIBOSOMAL_RNA (for --mode=ribo)"
      exit 1
    fi
  fi
fi

#	Also make sure $argGenome(.fa/.fasta) exists with a corresponding .bwt file
if [[ ! -e $argGenome ]]
then
  if [[ -e $argGenome.fa ]]
  then
    argGenome=$argGenome.fa
    echo "Using $argGenome as primary fasta file of genomic sequence"
  elif [[ -e $argGenome.fasta ]]
  then
    argGenome=$argGenome.fasta
    echo "Using $argGenome as primary fasta file of genomic sequence"
  else
    echo "ERROR: --genome does not specify a valid fasta file"
    echo "MISSING: $argGenome(.fa|fasta)"
    exit 1
  fi
fi

if [[ ! -e $argGenome.bwt ]]
then
  echo "ERROR: --genome points to a valid fasta file, but no BWA index"
  echo "MISSING: $argGenome.bwt"
  echo "Run bwa index $argGenome before running this step"
  exit 1
fi

if [[ $argPaired ]]
then
  echo "Sample is paired end, processing both fastq files together."
else
  echo "Sample is single end."
fi

# --unmapped - if unset, attempt to use mode-specific default, otherwise leave blank
if [[ (! $argUnmapped) && ($argMode == "ribo") && $RIBO_FILTERED_PATH ]]
then
  argUnmapped=$RIBO_FILTERED_PATH/$argBatch
  echo "Inferring --unmapped=$argUnmapped based on RIBO_FILTERED_PATH in project_info.sh"
elif [[ (! $argUnmapped) && ($argMode == "microbe") && $MICROBE_FILTERED_PATH ]]
then
  argUnmapped=$MICROBE_FILTERED_PATH/$argBatch
  echo "Inferring --unmapped=$argUnmapped based on MICROBE_FILTERED_PATH in project_info.sh"
elif [[ (! $argUnmapped) && ($argMode == "repeat") && REPEAT_FILTERED_PATH ]]
then
  argUnmapped=$REPEAT_FILTERED_PATH/$argBatch
  echo "Inferring --unmapped=$argUnmapped based on REPEAT_FILTERED_PATH in project_info.sh"
fi

# TO DO: Could have a similar default variable for --mapped, but currently leaving them in the BWA_PATH

# --mapped and --unmapped shouldn't be the same unless they're both blank
if [[ $argMapped && $argUnmapped ]]
then
  if [[ $argMapped == $argUnmapped ]]
  then
    echo "ERROR: --mapped and --unmapped must refer to different directories if both are set"
    exit 1
  fi
fi

# --threads - If not specified, check for DEFAULT_THREAD_COUNT, or default to 1
# TO DO: CHECK SLURM ENVIR BEFORE DEFAULT 
if [[ ! $argThreads ]]
then
  if [[ $SLURM_JOB_CPUS_PER_NODE ]]
  then
    argThreads=$SLURM_JOB_CPUS_PER_NODE
    echo "Inferring --threads=$argThreads based on number of CPUs assigned to slurm job."
  elif [[ $DEFAULT_THREAD_COUNT ]]
  then
    argThreads=$DEFAULT_THREAD_COUNT
    echo "Inferring --threads=$argThreads based on DEFAULT_THREAD_COUNT in project_info.sh."
  else
    argThreads=1
    echo "Setting --threads=$argThreads by default."
  fi
fi



# Check for dependencies and output the versions, e.g.
dependencies=("bwa" "samtools" "gzip" "gunzip")
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
    if [[ $tool == "samtools" || $tool == "bcftools" || $tool == "bwa" ]]
    then
      $tool 2>&1 | grep '^Version:'
    elif [[ $tool == "htseq-count" ]]
    then
      $tool 2>&1 | tail -n 3
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


# -- CORE TASK -- #

echo ""
echo "Running BWA task ($argMode mode) on sample: $argSample"
echo ""

# Any file added to this list will be deleted at the end
rmfiles=""
  
# Single-end sample processing
if [[ ! $argPaired ]]
then
  # Define full paths to input/output files based on Sample ID
  fastqin=$argInput/${argSample}$argInSuffix
  bwaout=$argOutput/$argSample
  if [[ ! $argMapped ]]
  then
    argMapped=$bwaout
  fi
  mappedout=$argMapped/${argSample}_mapped.fastq.gz
  
  if [[ ! $argUnmapped ]]
  then
    argUnmapped=$bwaout
  fi
  unmappedout=$argUnmapped/${argSample}_filtered.fastq.gz
  
  # Make sure we're not stomping on any existing output  
  if [[ -e $bwaout/output.bam ]]
  then
    echo "ERROR: BWA output already exists in $bwaout"
    echo "Remove the previous output before re-running this script"
    exit 1
  else
    echo "Primary alignment output will be in: $bwaout/output.bam"
  fi
  
  if [[ -e $mappedout ]]
  then
    echo "ERROR: Mapped reads file already exists at $mappedout"
    echo "Remove the previous output before re-running this script"
    exit 1
  else
    echo "Mapped reads will be output to: $mappedout"
  fi
  
  if [[ -e $unmappedout ]]
  then
    echo "ERROR: Unmapped reads file already exists at $unmappedout"
    echo "Remove the previous output before re-running this script"
    exit 1
  else
    echo "Unmapped reads will be output to: $unmappedout"
  fi

  # Check for .gz version if primary input file missing
  if [[ ! -e $fastqin ]]
  then
    if [[ -e $fastqin.gz ]]
    then
      # UPDATE 4/12/16 LJE
      # BWA can just take .fastq.gz file as input and auto-detect that it is compressed!
      # No need to decompress to temp file here!
      fastqin=$fastqin.gz
	  # --- OLD CODE --- #
	  # Need to unzip the file for bowtie
	  # echo "Decompressing fastq input to temporary copy"
	  #
	  # newfastq=$bwaout/${argSample}$argInSuffix
	  #
	  # echo "RUNNING: gunzip -c $fastqin.gz > $newfastq"
	  #
	  # if [[ ! $argDryRun ]]
	  # then
	  #   echo ""
	  #   echo -ne "Started at:\t"
	  #   date
	  #   echo ""
	  #  
	  #   mkdir -p $bwaout
	  #   gunzip -c $fastqin.gz > $newfastq
	  #   gunzip_exit_status=$?
	  #
	  #   echo ""
	  #   echo -ne "Finished at:\t"
	  #   date
	  #   echo ""
	  #  
  	  #   if [[ "$gunzip_exit_status" -ne "0" ]]
  	  #   then
  	  # 	echo "ERROR: gunzip returned non-zero exit status = $gunzip_exit_status"
  	  #	    exit $gunzip_exit_status
  	  #   fi
  	  #
  	  #   if [[ ! -e $newfastq ]]
  	  #   then
  	  #	    echo "ERROR: Can't find $newfastq after gunzip step"
  	  #	    exit 1
  	  #   fi
  	  #
  	  #   echo "Extracted fastq line count:"
  	  #   wc -l $newfastq
  	  #   echo ""
  	  # fi
  	  # fastqin=$newfastq
  	  # rmfiles=$fastqin 
  	  # --- END OLD CODE ---      
    else
  	  echo "ERROR: Missing $fastqin(.gz)"
  	  exit 1
    fi
  fi
  
  # (SKIPPING THIS) Setup appropriate RG header line
  # rgheader="ID:$argSample PL:ILLUMINA LB:L$lane SM:$sampleName PU:$flowcell"
  
  # --- BWA Alignment step --- #
  echo ""
  echo "RUNNING: bwa mem -v 2 -t $argThreads $argGenome $fastqin | samtools view -bh -o $bwaout/output.bam -"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
    
    mkdir -p $bwaout
    bwa mem -v 2 -t $argThreads $argGenome $fastqin | samtools view -bh -o $bwaout/output.bam -
    bwa_exit_status=$?
    
    # NOTE: Not bothering to sort the bam file here - that will only slow things down
    # Bonus: BWA should output the reads in the same order as the input fastq, so the unmapped file should be in same order
  	
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$bwa_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: bwa returned non-zero exit status = $bwa_exit_status"
  	  exit $bwa_exit_status
    fi
  
    if [[ ! -e $bwaout/output.bam ]]
    then
  	  echo "ERROR: bwa failed to produce alignment file"
  	  echo "MISSING: $bwaout/output.sam"
  	exit 1
    fi
  fi
  
  # --- Extract mapped and unmapped reads --- #
  echo ""
  echo "Extracting mapped and unmapped reads to fastq files"
  echo ""
  echo "RUNNING: samtools view -F 4 -b $bwaout/output.bam | samtools bam2fq - | gzip -c > $mappedout"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
  
    mkdir -p $argMapped
    samtools view -F 4 -b $bwaout/output.bam | samtools bam2fq - | gzip -c > $mappedout
    samtools_exit_status=$?
  
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$samtools_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: samtools or gzip returned non-zero exit status = $samtools_exit_status"
  	  exit $samtools_exit_status
    fi
  
    if [[ ! -e $mappedout ]]
    then
  	  echo "ERROR: Can't find mapped read fastq file"
  	  echo "MISSING: $mappedout"
  	  exit 1
    fi
  fi
  
  echo "RUNNING: samtools view -f 4 -b $bwaout/output.bam | samtools bam2fq - | gzip -c > $unmappedout"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
  
    mkdir -p $argUnmapped
    samtools view -f 4 -b $bwaout/output.bam | samtools bam2fq - | gzip -c > $unmappedout
    samtools_exit_status=$?
  
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$samtools_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: samtools or gzip returned non-zero exit status = $samtools_exit_status"
  	  exit $samtools_exit_status
    fi
  
    if [[ ! -e $unmappedout ]]
    then
  	  echo "ERROR: Can't find mapped read fastq file"
  	  echo "MISSING: $unmappedout"
  	  exit 1
    fi
  fi
else
  # Handle paired-end samples
  # Define full paths to input/output files based on Sample ID
  fastqinleft=$argInput/${argSample}_1$argInSuffix
  fastqinright=$argInput/${argSample}_2$argInSuffix
  bwaout=$argOutput/$argSample
  
  if [[ ! $argDryRun ]]
  then
    mkdir -p $bwaout
  fi
  
  if [[ ! $argMapped ]]
  then
    argMapped=$bwaout
  fi
  mappedoutleft=$argMapped/${argSample}_1_mapped.fastq.gz
  mappedoutright=$argMapped/${argSample}_2_mapped.fastq.gz
  
  if [[ ! $argUnmapped ]]
  then
    argUnmapped=$bwaout
  fi
  unmappedoutleft=$argUnmapped/${argSample}_1_filtered.fastq.gz
  unmappedoutright=$argUnmapped/${argSample}_2_filtered.fastq.gz
  
  # Make sure we're not stomping on any existing output  
  if [[ -e $bwaout/output.bam ]]
  then
    echo "ERROR: BWA output already exists in $bwaout"
    echo "Remove the previous output before re-running this script"
    exit 1
  else
    echo "Primary alignment output will be in: $bwaout/output.bam"
  fi
  
  if [[ (-e $mappedoutleft) || (-e $mappedoutright) ]]
  then
    echo "ERROR: Mapped reads file already exists at $mappedoutleft and/or $mappedoutright"
    echo "Remove the previous output before re-running this script"
    exit 1
  else
    echo "Mapped reads will be output to:"
    echo "$mappedoutleft"
    echo "$mappedoutright"
  fi
  
  if [[ (-e $unmappedoutleft) || (-e $unmappedoutright) ]]
  then
    echo "ERROR: Unmapped reads file already exists at $unmappedoutleft and/or $unmappedright"
    echo "Remove the previous output before re-running this script"
    exit 1
  else
    echo "Unmapped reads will be output to:"
    echo "$unmappedoutleft"
    echo "$unmappedoutright"
  fi
  
  # Check for .gz version if primary input files missing
  if [[ ! -e $fastqinleft ]]
  then
    if [[ -e $fastqinleft.gz ]]
    then
  	  # UPDATE 4/12/16 LJE
      # BWA can just take .fastq.gz file as input and auto-detect that it is compressed!
      # No need to decompress to temp file here!
      fastqinleft=$fastqinleft.gz
      # --- OLD CODE --- #
      # Need to unzip the file for bowtie
  	  # echo "Decompressing fastq input 1 to temporary copy"
      #  	
  	  # newfastq=$bwaout/${argSample}_1$argInSuffix
  	  #
  	  # echo "RUNNING: gunzip -c $fastqinleft.gz > $newfastq"
  	  #
  	  # if [[ ! $argDryRun ]]
  	  # then
  	  #   echo ""
  	  #   echo -ne "Started at:\t"
  	  #   date
  	  #   echo ""
  	  #
  	  #   gunzip -c $fastqinleft.gz > $newfastq
  	  #   gunzip_exit_status=$?
  	  #
  	  #   echo ""
  	  #   echo -ne "Finished at:\t"
  	  #   date
  	  #   echo ""
  	  #
  	  #   if [[ "$gunzip_exit_status" -ne "0" ]]
  	  #   then
  	  # 	echo "ERROR: gunzip returned non-zero exit status = $gunzip_exit_status"
  	  #	    exit $gunzip_exit_status
  	  #   fi
  	  #
  	  #   if [[ ! -e $newfastq ]]
  	  #   then
  	  #	    echo "ERROR: Can't find $newfastq after gunzip step"
  	  #	    exit 1
  	  #   fi
  	  #
  	  #   echo "Record of extracted fastq line count:"
  	  #   wc -l $newfastq
  	  #   echo ""
  	  # fi
  	  # fastqinleft=$newfastq
  	  # rmfiles=$fastqinleft   
  	  # --- END OLD CODE --- #   
    else
  	  echo "ERROR: Missing $fastqinleft(.gz)"
  	  exit 1
    fi
  fi
  
  if [[ ! -e $fastqinright ]]
  then
    if [[ -e $fastqinright.gz ]]
    then
      # UPDATE 4/12/16 LJE
      # BWA can just take .fastq.gz file as input and auto-detect that it is compressed!
      # No need to decompress to temp file here!
      fastqinright=$fastqinright.gz
      # --- OLD CODE --- #
  	  # Need to unzip the file for bowtie
  	  # echo "Decompressing fastq input 2 to temporary copy"
  	  # 
  	  # newfastq=$bwaout/${argSample}_2$argInSuffix
  	  # 
  	  # echo "RUNNING: gunzip -c $fastqinright.gz > $newfastq"
  	  # 
  	  # if [[ ! $argDryRun ]]
  	  # then
  	  #   echo ""
  	  #   echo -ne "Started at:\t"
  	  #   date
  	  #   echo ""
  	  # 
      #   gunzip -c $fastqinright.gz > $newfastq
      #   gunzip_exit_status=$?
  	  # 
  	  #   echo ""
  	  #   echo -ne "Finished at:\t"
  	  #   date
  	  #   echo ""
  	  # 
  	  #   if [[ "$gunzip_exit_status" -ne "0" ]]
  	  #   then
  	  # 	echo "ERROR: gunzip returned non-zero exit status = $gunzip_exit_status"
  	  # 	exit $gunzip_exit_status
  	  #   fi
  	  # 
  	  #   if [[ ! -e $newfastq ]]
  	  #   then
  	  # 	echo "ERROR: Can't find $newfastq after gunzip step"
  	  # 	exit 1
  	  #   fi
  	  # 
  	  #   echo "Record of extracted fastq line count:"
  	  #   wc -l $newfastq
  	  #   echo ""
  	  # fi
  	  # fastqinright=$newfastq
  	  # rmfiles="$rmfiles $fastqinright"
  	  # --- END OLD CODE --- #
    else
  	  echo "ERROR: Missing $fastqinright(.gz)"
  	  exit 1
    fi
  fi
  
  # (SKIPPING) Setup appropriate RG header line
  # rgheader="ID:$argSample PL:ILLUMINA LB:L$lane SM:$sampleName PU:$flowcell"
  
  # --- BWA Alignment step --- #
  echo ""
  echo "RUNNING: bwa mem -v 2 -t $argThreads $argGenome $fastqinleft $fastqinright | samtools view -bh -o $bwaout/output.bam -"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
    
    mkdir -p $bwaout
    bwa mem -v 2 -t $argThreads $argGenome $fastqinleft $fastqinright | samtools view -bh -o $bwaout/output.bam -
    bwa_exit_status=$?
    
    # NOTE: Not bothering to sort the bam file here - that will only slow things down
    # Bonus: BWA should output the reads in the same order as the input fastq, so the unmapped file should be in same order
  	
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$bwa_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: bwa returned non-zero exit status = $bwa_exit_status"
  	  exit $bwa_exit_status
    fi
  
    if [[ ! -e $bwaout/output.bam ]]
    then
  	  echo "ERROR: bwa failed to produce alignment file"
  	  echo "MISSING: $bwaout/output.sam"
  	  exit 1
    fi
  fi
  
  # --- Extract mapped and unmapped reads --- #
  echo ""
  echo "Extracting mapped and unmapped reads to fastq files"
  echo ""
  echo "RUNNING: samtools view -f 66 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $mappedoutleft"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
  
    mkdir -p $argMapped
    samtools view -f 66 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $mappedoutleft
    samtools_exit_status=$?
  
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$samtools_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: samtools or gzip returned non-zero exit status = $samtools_exit_status"
  	  exit $samtools_exit_status
    fi
  
    if [[ ! -e $mappedoutleft ]]
    then
  	  echo "ERROR: Can't find mapped read fastq file"
  	  echo "MISSING: $mappedoutleft"
  	  exit 1
    fi
  fi
  
  echo "RUNNING: samtools view -f 130 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $mappedoutright"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
  
    samtools view -f 130 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $mappedoutright
    samtools_exit_status=$?
  
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$samtools_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: samtools or gzip returned non-zero exit status = $samtools_exit_status"
  	  exit $samtools_exit_status
    fi
  
    if [[ ! -e $mappedoutright ]]
    then
  	  echo "ERROR: Can't find mapped read fastq file"
  	  echo "MISSING: $mappedoutright"
  	  exit 1
    fi
  fi
  
  echo "RUNNING: samtools view -F 2 -f 64 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $unmappedoutleft"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
  
    mkdir -p $argUnmapped
    samtools view -F 2 -f 64 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $unmappedoutleft
    samtools_exit_status=$?
  
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$samtools_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: samtools or gzip returned non-zero exit status = $samtools_exit_status"
  	  exit $samtools_exit_status
    fi
  
    if [[ ! -e $unmappedoutleft ]]
    then
  	  echo "ERROR: Can't find mapped read fastq file"
  	  echo "MISSING: $unmappedoutleft"
  	  exit 1
    fi
  fi
  
  echo "RUNNING: samtools view -F 2 -f 128 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $unmappedoutright"
  
  if [[ ! $argDryRun ]]
  then
    echo ""
    echo -ne "Started at:\t"
    date
    echo ""
  
    samtools view -F 2 -f 128 -u $bwaout/output.bam | samtools bam2fq - | gzip -c > $unmappedoutright
    samtools_exit_status=$?
  
    echo ""
    echo -ne "Finished at:\t"
    date
    echo ""
  
    if [[ "$samtools_exit_status" -ne "0" ]]
    then
  	  echo "ERROR: samtools or gzip returned non-zero exit status = $samtools_exit_status"
  	  exit $samtools_exit_status
    fi
  
    if [[ ! -e $unmappedoutright ]]
    then
  	  echo "ERROR: Can't find mapped read fastq file"
  	  echo "MISSING: $unmappedoutright"
  	  exit 1
    fi
  fi
fi

# --- Clean Up --- #
echo "Removing unmapped reads from $bwaout/output.bam"
echo "RUNNING: samtools view -F 4 -b $bwaout/output.bam > $bwaout/output.bam.aln ; mv $bwaout/output.bam.aln $bwaout/output.bam"
  
if [[ ! $argDryRun ]]
then
  echo ""
  echo -ne "Started at:\t"
  date
  echo ""
  
  samtools view -F 4 -b $bwaout/output.bam > $bwaout/output.bam.aln
  mv $bwaout/output.bam.aln $bwaout/output.bam
  
  echo ""
  echo -ne "Finished at:\t"
  date
  echo ""
fi

if [[ $rmfiles != "" ]]
then
  echo "Cleaning up temp files"
  echo ""
  echo "RUNNING: rm $rmfiles"

  if [[ ! $argDryRun ]]
  then
    rm $rmfiles
  fi
fi

echo "Finished BWA task."
echo ""
