#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 6/1/17
# 
# GOAL: Combine master sample tables for 2 data sets (treatment and control)
# Determine the lines that are complete in both sets, and drop excess/bad replicates
#
# Usage:
# Rscript combine_sample_tables.R [OPTIONS]
#  TRMTFILE= The table with sample information for treatment samples (default: sample_master_table.txt)
#  CTRLFILE= The table with sample information for control samples (default: ~/Resources/FlyBase/Dmel_r5.57_FB2014_03/DGRP/BaselineRNA/sample_master_table.txt)
#  TRMTCOL=  Column to add to designate treatment vs control samples (default: TRMT)
#  TRMTLABEL= Label for treatment samples (default: Treated)
#  CTRLLABEL= Label for control samples (default: Control)
#  OUTFILE=  Filename to output combined table (default: combined_master_table.txt)
