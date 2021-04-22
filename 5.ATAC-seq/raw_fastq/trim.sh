dir=/ifs1/User/shangyao/login123/ATAC-seq/fastq/clean_fastq/
cat fq.file| while read id
do
    echo $id
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    nohup trim_galore -q 26 --phred33 --length 40 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 &
done 
