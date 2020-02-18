# =========================================================================================
# = perform BLAT on NTR sequences and protein-coding sequences to see the characteristics =
# =========================================================================================

# 1. build ooc
# ============================================================

~/software/blat-v385/blat ~/flybase/fb-r5.57/dmel-all-chromosome-r5.57.fasta /dev/null /dev/null -makeOoc=bdgp5.ooc -repMatch=1024 > blat.ooc.log 2>&1 &

# 2. get rna sequences
# ============================================================

awk '$3 == "exon"' /home/ljeveret/Resources/FlyBase/Dmel_r5.57_FB2014_03/gff/dmel-all-transcriptome-r5.57-plus-r6.11-backport.gff > fb5.57.exon.gff
~/software/gffread-0.11.7.Linux_x86_64/gffread -w fb5.57.fasta -g ~/flybase/fb-r5.57/dmel-all-chromosome-r5.57.fasta fb5.57.exon.gff
~/software/gffread-0.11.7.Linux_x86_64/gffread -w ntr.fasta -g ~/flybase/fb-r5.57/dmel-all-chromosome-r5.57.fasta /home/ljeveret/Projects/DGRP_Baseline_RNAseq_Align/cuffmerge/merged_samples/cuffcompare.all.novel.gtf

cat fb5.57.fasta ntr.fasta > all.fasta

# 3. run blat
# ============================================================

~/software/blat-v385/blat ~/flybase/fb-r5.57/dmel-all-chromosome-r5.57.fasta all.fasta -ooc=bdgp5.ooc -out=blast8 all.blast8

# 4. get IDs
# ============================================================


