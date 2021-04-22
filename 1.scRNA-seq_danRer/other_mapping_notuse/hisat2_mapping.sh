cd clean_data
genome=/ifs1/User/shangyao/ncbi/public/sra/danRer_genome/hisat2_index/danRer
for id in `ls SRR*`
do
    prefix=$(echo $id|cut -d "_" -f 1)
    echo $prefix
    hisat2 -p 30  -x ${genome} -U $id -S ${prefix}.sam
    samtools view -S ${prefix}.sam -b > ${prefix}.bam ## 
    samtools sort ${prefix}.bam -o ${prefix}_sorted.bam ##by default sorted by chorosome position ,-n by reads name
    samtools index ${prefix}_sorted.bam ## produced th index files
    rm -f ${prefix}.sam ${prefix}.bam
done
