#!/bin/bash
set -e

#
# LJE - 6/14/16
#
# This script runs MMC on a single line means file
# The file should already be converted to CSV and filtered with filter_line_means.R
# Usage:
#  ./run_mmc.sh [Options] PATH/..._line_means.csv
# 
# Output will appear in PATH/ with similar file names to the input file
#
# Options:
#  -f --fast		- Automatically sets default resolution=51
#  -r --resolution=INT	- Sets the grid resolution for MMC algorithm (default 451, or 51 for fast)
#  -c --correlation=[pearson|spearman|kendall] - Correlation type used, default=pearson
#  -s --suffix=STRING	- If specified, this puts an extra string in the output file, AFTER default _mmc suffix for output files
#  --python=PATH 	- Path to desired version of python, default is my anaconda build of python: /home/ljeveret/Tools/anaconda2/bin/python
#  --mmc=PATH		- Full path to run-mmc.py script, default is my install: /home/ljeveret/Tools/pymmc-pipeline/run-mmc.py
#

# -- COMMAND LINE PARAMETERS -- #

# parse command line options
# NOTE: all variables specified here start with "arg"
args=`getopt -o "fr:c:s:" -l "fast,resolution:,correlation:,suffix:,python:,mmc:" -- "$@"`
echo "Running with command line arguments: $args"
eval set -- "$args"

while true;
do
  case $1 in
    -f|--fast)
      argFast=1
      shift 1;;
      
    -r|--resolution)
      argResolution=$2
      shift 2;;

    -c|--correlation)
      argCorrelation=$2
      shift 2;;
    
    -s|--suffix)
      argSuffix=$2
      shift 2;;
    
    --python)
      argPython=$2
      shift 2;;
    
    --mmc)
      argMMC=$2
      shift 2;;
    
    --)
      shift
      break;;
  esac
done

LINEMEANS=$1

# LINEMEANS param is required
if [[ "$LINEMEANS" == "" ]]
then
  echo "ERROR: Must specify a line means input file!"
  exit 1
fi

if [[ ! $argResolution ]]
then
  if [[ $argFast ]]
  then
    argResolution="51"
  else
    argResolution="451"
  fi
fi

if [[ ! $argCorrelation ]]
then
  argCorrelation="pearson"
fi

if [[ ("$argCorrelation" != "pearson") && ("$argCorrelation" != "spearman") && ("$argCorrelation" != "kendall") ]]
then
  echo "ERROR: --correlation=$argCorrelation is not valid, must be one of (pearson|spearman|kendall)"
  exit 1
fi

if [[ ! $argPython ]]
then
  argPython="/home/ljeveret/Tools/anaconda2/bin/python"
fi

if [[ ! $argMMC ]]
then
  argMMC="/home/ljeveret/Tools/pymmc-pipeline/run-mmc.py"
fi


# Check for dependencies and report the versions, e.g.
if [[ -e $argPython ]]
then
  echo "Found PYTHON: $argPython"
  $argPython --version
else
  echo "ERROR: Missing python executable at $argPython"
  exit 1
fi

if [[ -e $argMMC ]]
then
  echo "Found MMC: $argMMC"
else
  echo "ERROR: Missing MMC script at $argMMC"
  exit 1
fi


# -- MAIN TASK -- #

echo "Running MMC on $LINEMEANS"
FCOUNT=`wc -l $LINEMEANS | awk '{print $1-1}'`
echo "Running $argCorrelation MMC with grid-size=$argResolution on $FCOUNT features"
echo ""

# Remove the trailing .csv and replace "_line_means" ending with "_mmc"
OUTSTUB=`basename $LINEMEANS | sed 's/[.][^.]*$//'`
OUTSTUB=$(echo "$OUTSTUB" | sed 's/[_]line[_]means$//')"_mmc"
if [[ $argSuffix ]]
then
  OUTSTUB=$OUTSTUB"_"$argSuffix
fi

# Define corresponding output file names:
MMCOUT=$OUTSTUB".csv"
UNSMAP=$OUTSTUB"_unsorted_heatmap.png"
SRTMAP=$OUTSTUB"_sorted_heatmap.png"
SMTHMAP=$OUTSTUB"_smoothed_heatmap.png"

OUTPATH=`dirname $LINEMEANS`
echo "Output will be in $OUTPATH in files:"
echo "$MMCOUT"
echo "$UNSMAP"
echo "$SRTMAP"
echo "$SMTHMAP"
echo ""

# Now run the mmc script
# Several notes: The dependencies were fulfilled by setting up Anaconda, 
# so that is the default python executable that is used, but can be changed by command line
# Running MMC script on Hyperion also requires some minor configuration of matplotlib
# This is currently achieved by the matplotlibrc file in this directory,
# but it should probably be moved to the correct place in user directory to work in all project directories
echo "RUNNING:"
echo "$argPython $argMMC --verbose --correlation $argCorrelation --csv-in $LINEMEANS --sigma-num $argResolution --csv-out $MMCOUT --unsorted-heatmap $UNSMAP --sorted-heatmap $SRTMAP --smoothed-heatmap $SMTHMAP"

$argPython $argMMC --verbose --correlation $argCorrelation --csv-in $LINEMEANS --sigma-num $argResolution --csv-out $MMCOUT --unsorted-heatmap $UNSMAP --sorted-heatmap $SRTMAP --smoothed-heatmap $SMTHMAP

# TO DO: CHECK FOR ERRORS HERE?

echo ""
echo "Script completed successfully!"
