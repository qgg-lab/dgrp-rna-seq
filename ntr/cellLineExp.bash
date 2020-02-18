# kallisto index
~/qgg/software/kallisto-v0.46.1/kallisto index -i all.idx --make-unique all.fasta > kallisto.idx.log 2>&1 &

# quantification, download + kallisto
for srr in $(cat sra.acc)
do
	sbatch --export=srr="$srr" --output quant/"$srr".out --error quant/"$srr".err quant.sbatch
done

for srr in $(cat sra.acc)
do
  err=$(grep Error quant/"$srr".quant.log | wc -l)
  if [ $err -gt 0 ]
  then
  sbatch --export=srr="$srr" --output quant/"$srr".out --error quant/"$srr".err quant.sbatch
  fi
done
# summarize results
# gene ID
