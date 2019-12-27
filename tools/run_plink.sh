#!/bin/bash
#
# LJE - 2/23/17
# Very simple wrapper script to run plink
# Mainly created for easy submission to sbatch
#
# Usage:
#  run_plink.sh --genotypes=STRING --phenotypes=FILE --output=OUT [--prune]
#  The first three parameters are all currently REQUIRED
#  --genotypes is prefix to line genotype files, passed to --bfile param of plink,
#	these should first be created by running setup_plink.sh
#	e.g. freeze2.200line.common 
#  --phenotypes is the path to phenotype data file for mapping, passed to --pheno param of plink,
#   this file should be created by running filter_line_means.R with MODE=eqtl
#	e.g. expression/plink/combined_samples_gene_fpkm_VR_eQTL_F.pheno
#  --output is the path and prefix of all output files, passed to --out param of plink,
#	e.g. expression/plink/F
#  --prune (OPTIONAL) - adds an LD pruning step and changes the output to OUT.pruned
#     This is currently designed to be run IN PARALLEL with normal run, and can then be used
#	to filter SNPs AFTER the FDR procedure
#

# TO DO: Can have a default value for prefix based on --phenotypes (pull out the path and the last tag before .pheno)
# TO DO: Could try to autodetect the local freeze2.*.common files to get default for --genotypes
# TO DO: Could introduce other nice features that are standard in alignment pipeline scripts

# parse command line options
# NOTE: all variables specified here start with "arg"
args=`getopt -o "g:p:o:" -l "genotypes:,phenotypes:,output:,prune" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -g|--genotypes)
      argGenotypes=$2
      shift 2;;
    
    -p|--phenotypes)
      argPhenotypes=$2
      shift 2;;
    
    -o|--output)
      argOutput=$2
      shift 2;;
    
    --prune)
      argPrune=$1
      shift 1;;
    
    --)
      shift
      break;;
  esac
done

# Fill in default values
if [[ ! $argGenotypes ]]
then
  echo "ERROR: Missing required --genotypes parameter."
  exit 1
fi

if [[ ! $argPhenotypes ]]
then
  echo "ERROR: Missing required --phenotypes parameter."
  exit 1
fi

if [[ ! $argOutput ]]
then
  echo "ERROR: Missing required --output parameter."
  exit 1
fi

# Make sure genotype and phenotype files are present
if [[ ! -e $argGenotypes.bed ]]
then
  echo "ERROR: $argGenotypes.bed not found!"
  exit 1
fi

if [[ ! -e $argGenotypes.bim ]]
then
  echo "ERROR: $argGenotypes.bim not found!"
  exit 1
fi

if [[ ! -e $argGenotypes.fam ]]
then
  echo "ERROR: $argGenotypes.fam not found!"
  exit 1
fi

if [[ ! -e $argPhenotypes ]]
then
  echo "ERROR: $argPhenotypes not found!"
  exit 1
fi

# --prune option
# The actual parameters chosen for --indep-pairwise were suggested by Wen
# They made need to be further adjusted if it doesn't sufficiently remove the enrichment of SNPs within large inversions
if [[ $argPrune ]]
then
  argOutput=$argOutput.pruned
  echo "Activating LD pruning options, output files will start with $argOutput"
  pruneParams="--indep-pairwise 1000 10 0.5"
else
  pruneParams=""
fi

# Check what version of plink is to be used
echo "Using Plink installed at:"
which plink

echo ""


echo "Running:"
echo "plink --noweb --silent --bfile $argGenotypes --pheno $argPhenotypes --all-pheno --assoc --pfilter 0.00001 --out $argOutput $pruneParams"

plink --noweb --silent --bfile $argGenotypes --pheno $argPhenotypes --all-pheno --assoc --pfilter 0.00001 --out $argOutput $pruneParams

echo "Script completed successfully"
echo "See $argOutput.log for detailed Plink log."
