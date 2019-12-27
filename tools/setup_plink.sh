#!/bin/bash
set -e

#
# LJE - 2/23/17
#
# Steps to setup common plink files for all feature analysis
#
# Usage:
#  setup_plink.sh
# 
# Output will appear in PATH/ with similar file names to the input file
#
# Options:
#  -s --samples=FILE	- Path to sample_master_table.txt file
#  -g --genotypes=FILE	- Path to genotypes file (minus .bed/.fam suffix)
#  -p --prefix=STRING	- Prefix to put on front of all output files (default=freeze2)
#

# -- COMMAND LINE PARAMETERS -- #

# parse command line options
# NOTE: all variables specified here start with "arg"
args=`getopt -o "s:g:p:" -l "samples:,genotypes:,prefix:" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -s|--samples)
      argSamples=$2
      shift 2;;

    -g|--genotypes)
      argGenotypes=$2
      shift 2;;
    
    -p|--prefix)
      argPrefix=$2
      shift 2;;
    
    --)
      shift
      break;;
  esac
done

# Fill in default values
if [[ ! $argSamples ]]
then
  argSamples="sample_master_table.txt"
  echo "Setting --samples=$argSamples by default."
fi

if [[ ! $argGenotypes ]]
then
  argGenotypes="/home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/freeze2.jgil.hpp.biallelic"
  echo "Setting --genotypes=$argGenotypes by default."
fi

if [[ ! $argPrefix ]]
then
  argPrefix="freeze2"
  echo "Setting --prefix=$argPrefix by default."
fi

# Make sure $argSamples and $argGenotypes exist
if [[ ! -e "$argSamples" ]]
then
  echo "ERROR: $argSamples not found!"
  exit 1
fi

GENOBEDFILE="$argGenotypes.bed"
if [[ ! -e $GENOBEDFILE ]]
then
  echo "ERROR: $GENOBEDFILE not found!"
  exit 1
fi


# Check what version of plink is to be used
echo "Using Plink installed at:"
which plink

echo ""


# -- MAIN TASK -- #

# Create plink ID file with row for each DGRP line in the format:
# line_XXX line_XXX
# line_YYY line_YYY
# ...
# (Space separated)
LNUM=`cut -f 5 $argSamples | grep -v 'LINE' | sort | uniq | wc -l`
PREFIX="$argPrefix.${LNUM}line"
PLINKFILE="$PREFIX.id.plink"
echo "$argSamples contains samples from $LNUM lines, storing line ID list in $PLINKFILE"
cut -f 5 $argSamples | grep -v 'LINE' | sort | uniq | awk '{print "line_"$1" line_"$1}' > $PLINKFILE
echo ""

# extract common genotypes
PLINKOUT="$argPrefix.${LNUM}line.common"
echo "Extracting genotypes for common alleles in $LNUM lines from $argGenotypes"
plink --noweb --silent --bfile $argGenotypes --keep $PLINKFILE --geno 0.25 --maf 0.05 --make-bed --out $PLINKOUT
echo "Extraction complete - see $PREFIX.common.log for Plink version, processing details."
echo ""

SNPFILE="$PLINKOUT.snp"
cut -f 2 $PLINKOUT.bim > $SNPFILE
SNPCOUNT=`wc -l $SNPFILE | awk '{print $1}'`

echo "There are $SNPCOUNT SNPs that can be analyzed as eQTLs in your $LNUM lines, see list of SNP IDs in $SNPFILE"

echo "Script completed successfully!"
