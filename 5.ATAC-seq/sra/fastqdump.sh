for i in *.sra;
do
    echo $i
    fastq-dump --gzip --split-3 $i
done
