dir=~/ncbi/public/sra/fastq
cd $dir
outdir=~/ncbi/public/sra/clean_data
ls *.gz | while read id
do
    echo $id
    nohup trim_galore -q 26 --phred33 --length 10 --gzip -e 0.1 --stringency 3  -o ${outdir} ${dir}/${id} &
done 
