for id in `ls *.sra`
do
    nohup fastq-dump --split-files --gzip $id >$id.log 2>&1 &
done
