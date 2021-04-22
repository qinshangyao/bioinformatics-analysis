for i in *.sra
do
    echo running $i
    fastq-dump --gzip $i
done
