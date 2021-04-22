tail -11 config.clean > config.clean2
##should be obsolute dictionary
bowtie2_index=/ifs1/User/shangyao/genome/bowtie2_index/mm10
## 一定要搞清楚自己的bowtie2软件安装在哪里，以及自己的索引文件在什么地方！！！
#source activate atac 
cat config.clean2 |while read id;
do 
	echo $id
	arr=($id)
	fq2=${arr[2]}
	fq1=${arr[1]}
	sample=${arr[0]}
	#### mapping start!!
	bowtie2  -p 32  --very-sensitive -X 2000 -x  ${bowtie2_index} -1 $fq1 -2 $fq2 |samtools sort  -O bam  -@ 32 -o - > ${sample}.raw.bam 
	samtools index ${sample}.raw.bam 
	bedtools bamtobed -i ${sample}.raw.bam  > ${sample}.raw.bed
	samtools flagstat ${sample}.raw.bam  > ${sample}.raw.stat
done
