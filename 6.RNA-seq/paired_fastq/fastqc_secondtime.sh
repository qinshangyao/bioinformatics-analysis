for id in `ls *fastq`
do
    fastqc $id -o ./fastqc_result -q
done
