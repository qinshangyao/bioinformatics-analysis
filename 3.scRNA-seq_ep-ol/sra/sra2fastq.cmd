for id in `ls SRR1151094*.sra`
do
    nohup fastq-dump --split-files --gzip $id >$id.log 2>&1 &
done
