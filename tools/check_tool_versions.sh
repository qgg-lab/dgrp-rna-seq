#!/usr/bin/bash

echo "SAMTools (should be v1.1):"
which samtools
samtools --version
echo ""

echo "CutAdapt (should be v1.6):"
which cutadapt
cutadapt --version
echo ""

echo "BWA (should v0.7.*):"
which bwa
bwa 2>&1 | head -n 4
echo ""
