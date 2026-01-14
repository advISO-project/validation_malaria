indir=$1
outfile=fastq_readcounts.csv
echo "file,count">${outfile} 
for f in ${indir}/*_1.fastq.gz; do
    fr=`realpath ${f}`
    base=`basename ${fr}`
    count=$(($(zcat ${fr}| wc -l) / 4))
    echo "${base},${count}"
done >> ${outfile}  