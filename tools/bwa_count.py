#!/usr/bin/python
#
# LJE - 12/15/15
#
# GOAL: Given a SAM/BAM file from BWA, count how many reads aligned to each target sequence.
# For reads with multiple alignments, count the read as 1/N for each alignment
# This is primarily intended to quantify alignments to entries in contaminant databases (rRNA, microbes)
# If multiple bam files are given, each one is read, and total counts across all bam files are given
#
# Usage:
# ./bwa_count.py [Options] align1.bam [align2.bam ...] > counts.txt
#
# DEPENDS ON: pysam
# NOTE: Using API for v0.8.0 b/c that's the version installed under default Python environment on hyperion
# For notes on what has changed in the API since then, see: http://pysam.readthedocs.org/en/latest/release.html#release-0-8-1
#
# UPDATE 2/9/16 - Added --nosplit option to count reads that align to distinct sets of target sequences
#
# UPDATE 10/18/16 - Added --mask option to provide a list of reads to mask
#
# UPDATE 10/26/16 - Added --rmdup option to exclude reads marked as duplicates
#
# TO DO: Drop additional alignments for a read if the difference in MAQ scores is above a certain threshold
#

# -- REQUIRED MODULES -- #

# Import into main namespace:
from sys import *
from optparse import OptionParser
# from datetime import datetime

# Import in their own namespace:
import os.path
import os
import pysam
import fileinput
# import re


appname = os.path.basename(argv[0])
appcall = " ".join(argv)

# -- COMMAND LINE PARAMETERS -- #

parser = OptionParser()

parser.add_option("--nosplit", dest="nosplit", action="store_true",
	help="For reads that align to multiple targets, report the set of target sequences instead of splitting the count")
parser.add_option("-m", "--mask", dest="maskReadFile",
                  help="FILE containing list of read IDs to mask", metavar="FILE")
parser.add_option("--rmdup", dest="rmdup", action="store_true",
	help="Skip reads flagged as duplicates when counting")
parser.add_option("--verbose", dest="verbose", action="store_true",
	help="Verbose progress reporting and warnings")

(options, args) = parser.parse_args()

# Remaining args are assumed to be list of bam files to parse


# -- GLOBAL DATA STRUCTURES -- #

# This dictionary will be keyed on read ID
# ONLY read IDs to mask should go in here at all
# Each value is an integer indicating the number of times the read ID was actually masked
# so 0 = read ID was in mask list but was never seen in primary input file
maskReadDict = {}

# This dictionary will be keyed on read ID
# Each value is a list of the observed target sequences from each alignment for that read
readTargetDict = {}

# This dictionary will be keyed on each target sequence
# Each value will be the weighted count of aligned reads
targetCountDict = {}

# Count number of lines handled
totalSamReads = 0
totalDupReads = 0
totalFiles = 0


# -- SUBROUTINES -- #

# Open and loop through a single text file that lists read IDs to mask (one per line)
def handleMaskFile(maskFile):
	global maskReadDict
	# Make sure the file exists
	if os.path.isfile(maskFile):
		# Open the file, add each read ID to maskReadDict
		for line in fileinput.input(maskFile):
			# Strip out surrounding white space and line breaks
			read = line.strip(" \t\n")
			# If not blank, add to maskReadDict
			if read != "" and read not in maskReadDict:
				maskReadDict[read] = 0
	else:
		# Return -1 error code
		return -1


# Add a target (reference) name to readTargetDict under key = read
# Writes directly to readTargetDict global dictionary
def addReadTarget(read, target):
	global readTargetDict
	if read in readTargetDict:
		# Append target to existing list
		readTargetDict[read].append(target)
	else:
		# Start new list
		readTargetDict[read] = [target]


# Handle a single line of data from SAM format input
# Also has to pass the alignment file object in order to call getrname() to look up target sequence name by index
# samRead should be of type pysam.AlignedSegment
# Return Value: 1 if successful, 0 if line empty or can't parse
def handleSamLine(samRead, parentFile):
	global maskReadDict
	global totalDupReads
	global options
	# Check if read has an alignment
	if samRead.rname != -1:
		# Yes, now check if readName is in maskReadDict
		readName = samRead.qname
		if readName in maskReadDict:
			# Masked, just increment masking count for this read
			maskReadDict[readName] += 1
		else:
			# Not masked, check if it's a duplicate when rmdup = True
			if options.rmdup and samRead.is_duplicate:
				# Read is duplicate, increment the counter
				totalDupReads += 1
			else:
				# Not a duplicate, add read,target pair to readTargetDict
				targetName = parentFile.getrname(samRead.rname)
				addReadTarget(readName, targetName)
				# Loop over tags to check for 'XA'
				for srTag in samRead.tags:
					if srTag[0] == 'XA':
						readXAstr = srTag[1]
						# First split on ';' to loop over multiple entries
						alignments = readXAstr.split(';')
						for aln in alignments:
							# Ignore empty strings
							if aln != '':
								# Split the XA field on ',', first entry is the target sequence ID
								alnTarget = aln.split(',')[0]
								# Again, ignore if empty
								if alnTarget != '':
									# Add to readTargetDict
									addReadTarget(readName, alnTarget)
		return 1
	else:
		return 0


# Open and loop through a single SAM/BAM file
# Return Value: Number of lines handled (or -1 if there's a problem)
def handleSamFile(samfile):
	# Make sure the file exists before starting sam
	if os.path.isfile(samfile):
		# Open samfile with pysam
		# (Using deprecated API "Samfile" instead of "AlignmentFile" for backwards compatibility)
		samin = pysam.Samfile(samfile)
		# Iterate over all aligned reads, until_eof=True avoids errors when there is no index
		readCount = 0
		for sRead in samin.fetch(until_eof=True):
			readCount += handleSamLine(sRead, samin)
		samin.close()
		return readCount
	else:
		# Return -1 error code
		return -1


# -- MAIN CODE -- #

# Check for mask file
if options.maskReadFile is not None:
	handleMaskFile(options.maskReadFile)
	if options.verbose:
		print >> stderr, appname + ": Read " + str(len(maskReadDict)) + " masked read IDs from " + options.maskReadFile

# Loop over each SAM/BAM file in args
for fname in args:
	if options.verbose:
		print >> stderr, appname + ": Reading " + fname
	hReturn = handleSamFile(fname)
	# Check return value
	if hReturn > 0:
		# If return value is positive integer, it's the count of lines parsed
		if options.verbose:
			print >> stderr, appname + ": Parsed " + str(hReturn) + " lines from " + fname
		totalSamReads += hReturn
		totalFiles += 1
	elif hReturn == 0:
		if options.verbose:
			print >> stderr, appname + ": " + fname + " was empty or could not be parsed."
	else:
		# Error code, report this even in non-verbose mode
		print >> stderr, appname + " WARNING: There was a problem reading " + fname

# Generate target->count dictionary from read->target list dictionary
# Loop over every key in rtDict
for read in readTargetDict:
	if options.nosplit:
		# In --nosplit mode, the key is a comma-separated sorted list of all unique targets
		target = ",".join(sorted(set(readTargetDict[read])))
		if target in targetCountDict:
			targetCountDict[target] += 1
		else:
			targetCountDict[target] = 1
	else:
		# Default mode: Each read counts as '1' for all alignments
		# So weight this read as 1/N if there are multiple alignments
		readWt = 1.0 / len(readTargetDict[read])
		# Loop over targets here
		for target in readTargetDict[read]:
			if target in targetCountDict:
				targetCountDict[target] += readWt
			else:
				targetCountDict[target] = readWt

# Output contents of targetCountDict to stdout
for target in targetCountDict:
	print >> stdout, target + "\t" + str(targetCountDict[target])

# (Verbose) Report number of files and alignment lines handled, total number of read and target IDs seen
if options.verbose:
	print >> stderr, appname + " COMPLETE: Read " + str(totalSamReads) + " alignment entries from " + str(totalFiles) + " SAM/BAM file(s)"
	maskSeenCnt=0
	maskMissCnt=0
	for read in maskReadDict:
		if maskReadDict[read] == 0:
			maskMissCnt += 1
		else:
			maskSeenCnt += 1
	if maskSeenCnt > 0:
		print >> stderr, appname + " COMPLETE: Masked " + str(maskSeenCnt) + " reads from being counted"
	if maskMissCnt > 0:
		print >> stderr, appname + " NOTE: " + str(maskMissCnt) + " reads were in mask list but not in SAM file(s)"
	if totalDupReads > 0:
		print >> stderr, appname + " COMPLETE: Did not count " + str(totalDupReads) + " reads that were flagged as duplicates"
	print >> stderr, appname + " COMPLETE: Processed " + str(len(readTargetDict)) + " reads aligning to " + str(len(targetCountDict)) + " target sequences"
	print >> stderr, appname + " NOTE: Input files may contain multiple alignments per entry, and may contain multiple entries for the same read"

# Done!
