#!/bin/bash

# LJE - 3/25/15
#
# run_star.sh - Script for running STAR on sequence files for DGRP Baseline RNA-seq project
# This does NOT handle submitting batch jobs to slurm controller - use in conjunction with run_batch_task.sh for that
#
# Usage:
# ./run_star.sh [options] --batch=batchname --sample=SampleID
# 
# Options:
#  --sample=sampleID - Refers to the main part of input file name (minus certain standard suffix portions, see below), also used to name the output files
#		Ultimately, the full input path is: inputpath/(batch)/sampleID(suffix).fastq(.gz)
#  --input=path - Path to the input directory, default is to first search for environment variables in this order:
#		$REPEAT_FILTERED_PATH, $MICROBE_FILTERED_PATH, $RIBO_FILTERED_PATH, $TRIMMED_FASTQ_PATH, or $PROC_FASTQ_PATH
#		If none are defined in project_info.pl, defaults to ./
#  --suffix=string - The "_filtered" or "_trimmed" portion of input fastq(.gz) files. Default depends on inferred input path, can be blank
#  --output=path - Path to the output directory, default is $STAR_PATH, as defined in project_info.sh (if missing, checks ALIGNED_PATH, then defaults to STAR/)
#		The actual output goes in a subdirectory named sampleID in the path indicated here
#		If STAR output for this sample exists already, this script will quit with an error
#  --batch=batchname - Convenience option, adds a matching subdirectory to both the input and output paths
#  --genome=path - Path to STAR genome index directory. Default checks $FLY_GENOME_STAR, then $FLY_GENOME, then hits an error...
#  --mismatch=NUMERIC - Parameter to adjust --outFilterMismatchNoverLmax (default = 0.04)
#     NOTE: This is the rate of mismatches per length of read, so for this allows 4 mismatches for a 100bp read, 5 for a 125bp read
#  --paired - Set this option when sample contains paired end reads
#  --unstranded - Adds options appropriate for unstranded data (mainly "--outSAMstrandField intronMotif", so the BAM file can be used with cufflinks)
#  --lane - For BAM file header only, skipped if blank
#  --name - source sample name for BAM file header, skipped if blank
#  --flowcell = Flow cell ID for BAM file header, skipped if blank
#  --threads=N - Specifies the number of threads to use, if unspecified, will check SLURM environment, then DEFAULT_THREAD_COUNT, otherwise defaults to 1
#  --dryrun - Parameter indicating that task should not actually be run, but any pre-task checks should be run,
#		and the final formatted command should be echo'd. This is useful for testing and debugging.
#		Basically, if dryrun is specified, there should be no output or file manipulation (code below must honor this though)
#
#

# parse command line options
args=`getopt -o "s:i:o:b:g:m:t:pud" -l "sample:,input:,suffix:,output:,batch:,genome:,mismatch:,paired,unstranded,lane:,name:,flowcell:,threads:,dryrun" -- "$@"`
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
    
    --suffix)
      argSuffix=$2
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
    
    -m|--mismatch)
      argMismatch=$2
      shift 2;;
    
    -p|--paired)
      argPaired=1
      shift 1;;
    
    -u|--unstranded)
      argUnstranded=1
      shift 1;;
    
    --lane)
      argLane=$2
      shift 2;;
    
    --name)
      argName=$2
      shift 2;;
    
    --flowcell)
      argFlowcell=$2
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
  echo "No project_info.sh, attempting to proceed without ."
fi

# Check each parameter
# Specify required params here and set defaults for others

# --sample - Required
if [[ ! $argSample ]]
then
  echo "ERROR: Must specify --sample parameter"
  exit 1
fi

# --input, if not specified check for $MICROBE_FILTERED_PATH, $RIBO_FILTERED_PATH, $TRIMMED_FASTQ_PATH, or $PROC_FASTQ_PATH, or default to ./
# Also sets corresponding argSuffix here if not specified on command line
if [[ ! $argInput ]]
then
  if [[ $REPEAT_FILTERED_PATH ]]
  then
    argInput=$REPEAT_FILTERED_PATH
    echo "Inferring --input=$argInput based on REPEAT_FILTERED_PATH in project_info.sh"
    if [[ ! $argSuffix ]]
    then
      argSuffix="_filtered"
      echo "Using --suffix=$argSuffix by default"
    fi
  elif [[ $MICROBE_FILTERED_PATH ]]
  then
    argInput=$MICROBE_FILTERED_PATH
    echo "Inferring --input=$argInput based on MICROBE_FILTERED_PATH in project_info.sh"
    if [[ ! $argSuffix ]]
    then
      argSuffix="_filtered"
      echo "Using --suffix=$argSuffix by default"
    fi
  elif [[ $RIBO_FILTERED_PATH ]]
  then
    argInput=$RIBO_FILTERED_PATH
    echo "Inferring --input=$argInput based on RIBO_FILTERED_PATH in project_info.sh"
    if [[ ! $argSuffix ]]
    then
      argSuffix="_filtered"
      echo "Using --suffix=$argSuffix by default"
    fi
  elif [[ $TRIMMED_FASTQ_PATH ]]
  then
    argInput=$TRIMMED_FASTQ_PATH
    echo "Inferring --input=$argInput based on TRIMMED_FASTQ_PATH in project_info.sh"
    if [[ ! $argSuffix ]]
    then
      argSuffix="_trimmed"
      echo "Using --suffix=$argSuffix by default"
    fi
  elif [[ $PROC_FASTQ_PATH ]] 
  then
    argInput=$PROC_FASTQ_PATH
    echo "Inferring --input=$argInput based on PROC_FASTQ_PATH in project_info.sh"
    if [[ ! $argSuffix ]]
    then
      argSuffix="_trimmed"
      echo "Using --suffix=$argSuffix by default"
    fi
  else
    argInput="."
    echo "Setting --input=$argInput by default"
    if [[ ! $argSuffix ]]
    then
      argSuffix=""
      echo "Leaving --suffix blank by default."
    fi
  fi
fi

# --output, if not specified check for STAR_PATH or ALIGNED_PATH, or default to STAR/
if [[ ! $argOutput ]]
then
  if [[ $STAR_PATH ]]
  then
    argOutput=$STAR_PATH
    echo "Inferring --output=$argOutput based on STAR_PATH in project_info.sh"
  elif [[ $ALIGNED_PATH ]]
  then
    argOutput=$ALIGNED_PATH
    echo "Inferring --output=$argOutput based on ALIGNED_PATH in project_info.sh"
  else
    argOutput="STAR"
    echo "Setting --output=$argInput by default"
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

# --genome=path, if not specified check $FLY_GENOME_STAR, then $FLY_GENOME (one of these required)
#	Also make sure $argGenome/SA exists
if [[ ! $argGenome ]]
then
  if [[ $FLY_GENOME_STAR ]]
  then
    argGenome=$FLY_GENOME_STAR
    echo "Inferring --genome=$argGenome from FLY_GENOME_STAR in project_info.sh"
  elif [[ $FLY_GENOME ]]
  then
    argGenome=$FLY_GENOME
    echo "Inferring --genome=$argGenome from FLY_GENOME in project_info.sh"
  else
    echo "ERROR: Must specify genome index with --genome, or in project_info.sh with FLY_GENOME_STAR"
    exit 1
  fi
fi

if [[ ! -e $argGenome/SA ]]
then
  echo "ERROR: --genome does not specify a valid STAR genome index"
  echo "MISSING: $argGenome/SA"
  exit 1
fi

# --mismatch - Default to 0.04 if not set
if [[ ! $argMismatch ]]
then
  argMismatch="0.04"
  echo "Setting --mismatch=$argMismatch by default."
fi

# --paired - Just report if this parameter has been set
if [[ $argPaired ]]
then
  echo "Processing paired-end read libraries"
fi

# --unstranded - Just report if this parameter has been set
if [[ $argUnstranded ]]
then
  echo "Processing unstranded read libraries"
fi

# --threads - If not specified, check for DEFAULT_THREAD_COUNT, or default to 1
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

# --dryrun - Just report if this is the case
if [[ $argDryRun ]]
then
  echo "---- DRY RUN ----"
  echo "Tasks will not actually run, there will be no output!"
fi


# Check for dependencies and output the versions, e.g.
dependencies=("zcat" "STAR")
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

if [[ ! -d $argInput ]]
then
  echo "ERROR: Input directory does not exist"
  echo "MISSING: $argInput"
  exit 1
fi


# -- CORE TASK -- #

# Check for input files (differs for paired and single reads)
if [[ $argPaired ]]
then
  fastqinleft=$argInput/${argSample}_1${argSuffix}.fastq
  fastqinright=$argInput/${argSample}_2${argSuffix}.fastq
  # Check for .gz version if primary input file missing
  # Because of the --readFilesCommand param, if the first read file is gzipped, the other one must be as well
  if [[ ! -e $fastqinleft ]]
  then
    if [[ (-e $fastqinleft.gz) && (-e $fastqinright.gz) ]] 
    then
      # The additional --readFilesCommand zcat param is required to tell STAR how to decompress the gzipped files to stdout
      fastqin="$fastqinleft.gz $fastqinright.gz --readFilesCommand zcat"
    else
      echo "ERROR: Missing either $fastqinleft(.gz) or $fastqinright(.gz)"
      echo "Note that if one input file is compressed, the corresponding mate pair file must also be compressed."
      exit 1
    fi
  else
    # Input files are not zipped
    fastqin="$fastqinleft $fastqinright"
  fi
  # Regardless of whether zipped or not, add the --alignMatesGapMax param for paired reads
  fastqin="$fastqin --alignMatesGapMax 1000000"
else
  fastqin=$argInput/${argSample}${argSuffix}.fastq
  # Check for .gz version if primary input file missing
  if [[ ! -e $fastqin ]]
  then
    if [[ -e $fastqin.gz ]]
    then
      # The additional --readFilesCommand zcat param is required to tell STAR how to decompress the gzipped file to stdout
      fastqin="$fastqin.gz --readFilesCommand zcat"
    else
      echo "ERROR: Missing $fastqin(.gz)"
      exit 1
    fi
  fi
fi
# At the end of this, the key variable is still fastqin regardless of paired/single
# This variable contains all parameters relevant to specifying the input file(s) and related params
  
  
# Check for pre-existing output
starout=$argOutput/$argSample
if [[ -e $starout/Aligned.sortedByCoord.out.bam ]]
then
  echo "ERROR: STAR output already exists in $starout"
  echo "Remove the previous output before re-running this script"
  exit 1
fi
  
if [[ -e $starout/align_summary.out ]]
then
  echo "ERROR: Intended output directory already contains Tophat2 alignment results"
  echo "Changed ALIGNED_PATH in project_info.sh before proceeding"
  exit 1
fi

# U R HERE - 3/25/15 - Need to add parameters for these header components - allow for default values
  
# Setup appropriate RG header line, including optional fields
rgheader="ID:$argSample PL:ILLUMINA"
if [[ $argLane ]]
then
  rgheader="$rgheader LB:L$argLane"
fi
if [[ $argName ]]
then
  rgheader="$rgheader SM:$argName"
fi
if [[ $argFlowcell ]]
then
  rgheader="$rgheader PU:$argFlowcell"
fi
 
# Build out the string of STAR parameters
starParams="--runThreadN $argThreads --outFilterMismatchNoverLmax $argMismatch --outFilterIntronMotifs RemoveNoncanonicalUnannotated"
if [[ $argUnstranded ]]
then
  starParams="$starParams --outSAMstrandField intronMotif"
fi
starParams="$starParams --outSAMattrRGline $rgheader --genomeDir $argGenome --readFilesIn $fastqin --outFileNamePrefix $starout/ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --limitBAMsortRAM 10000000000"
   
echo ""
echo "RUNNING: STAR $starParams"
  
if [[ ! $argDryRun ]]
then
  echo ""
  echo -ne "Started at:\t"
  date
  echo ""
    
  mkdir -p $starout
  STAR $starParams
  star_exit_status=$?
  # TO DO: Look over 2-pass mode documentation - best way to run this is to collect all novel splice junctions, filter to high-conf ones, rebuild genome index, and run all samples again
  # Could presumably collect all high-confidence novel-transcripts at this point as well
  # NOTE: There is also the --outFilterMismatchNmax parameter that directly sets the number of mismatches allowed (analogous to tophat option), 
  #	but I think scaling it to the aligned portion of the read makes more sense.  The setting above will allow at most 4 mismatches for 100bp of aligned sequence
  # NOTE: The documentation recommended this parameter: "--outFilterIntronMotifs RemoveNoncanonical" when using Cufflinks,
  #	If I understand things correctly, "non-canonical" splice junctions are those that don't conform to typical splicing motifs
  #   I changed the option to only filter out splice junctions with that are both non-canonical AND unannotated (if they're in the FlyBase transcriptome I don't see a reason to filter them out)
  # NOTE: Found the option to output unmapped reads back to fastq format: --outReadsUnmapped Fastx
  #	This is definitely useful for doing rRNA filtering, but also useful in general to do assembly/microbe checks where there is large number of unmapped reads
  #	It might eventually be useful to have an extra step at the end to move this file somewhere else and/or compress it, or delete it
  # NOTE: Also has options to process the BAM file into Wig/BedGraph output tracks (requires a separate invokation of STAR in different mode)
  # NOTE: There are a number of --outSJ* parameters that will be useful if building novel SJDB and index for 2-pass approach. 
  #	Filtering can also be done on the file itself as it contains info about each splice junction in tab-delimited format
  # NOTE: Looks to me like STAR is using 0.9Gb per thread for Dmel genome, doesn't need nearly as much memory as for human genome
  # NOTE: In terms of timing, this only took 5 minutes with 4 threads,
  #	The same fastq file, same number of threads took Tophat ~33m, so that's more than 6-fold speed up...
  # NOTE: For some reason doing the rRNA alignment is MUCH slower, I may need to increase the SA size when building the index
  # 	Also, the amount of memory allowed for BAM sorting is based on genome size, so the --limitBAMsortRAM parameter helps ensure there's enough memory during the rRNA alignment. May need to increase further for bigger samples
  #
	    
  echo ""
  echo -ne "Finished at:\t"
  date
  echo ""
    
  if [[ "$star_exit_status" -ne "0" ]]
  then
    echo "ERROR: STAR returned non-zero exit status = $star_exit_status"
    exit $star_exit_status
  fi
    
  if [[ ! -e $starout/Aligned.sortedByCoord.out.bam ]]
  then
    echo "ERROR: STAR failed to produce alignment file"
    echo "MISSING: $starout/Aligned.sortedByCoord.out.bam"
    exit 1
  fi
fi
