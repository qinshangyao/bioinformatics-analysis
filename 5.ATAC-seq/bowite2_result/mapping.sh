ls ~/login123/ATAC-seq/clean_fastq/*_1.fq.gz >1
ls ~/login123/ATAC-seq/clean_fastq/*_2.fq.gz >2
ls ~/login123/ATAC-seq/clean_fastq/*_2.fq.gz |cut -d "/" -f 8 | cut -d "_" -f 1 >0
paste 0 1 2 > config.clean
rm -f 0 1 2
##should be obsolute dictionary
bowtie2_index=/ifs1/User/shangyao/genome/bowtie2_index/mm10
## 一定要搞清楚自己的bowtie2软件安装在哪里，以及自己的索引文件在什么地方！！！
#source activate atac 
cat config.clean |while read id;
do 
	echo $id
	arr=($id)
	fq2=${arr[2]}
	fq1=${arr[1]}
	sample=${arr[0]}
	#### mapping start!!
	bowtie2  -p 10  --very-sensitive -X 2000 -x  ${bowtie2_index} -1 $fq1 -2 $fq2 |samtools sort  -O bam  -@ 10 -o - > ${sample}.raw.bam 
	samtools index ${sample}.raw.bam 
	bedtools bamtobed -i ${sample}.raw.bam  > ${sample}.raw.bed
	samtools flagstat ${sample}.raw.bam  > ${sample}.raw.stat
	
	# https://github.com/biod/sambamba/issues/177
	sambamba markdup --overflow-list-size 600000  --tmpdir='./'  -r ${sample}.raw.bam  ${sample}.rmdup.bam
	samtools index   ${sample}.rmdup.bam 
	
	## ref:https://www.biostars.org/p/170294/ 
	## Calculate %mtDNA:
	mtReads=$(samtools idxstats  ${sample}.rmdup.bam | grep 'chrM' | cut -f 3)
	totalReads=$(samtools idxstats  ${sample}.rmdup.bam | awk '{SUM += $3} END {print SUM}')
	echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'
	
	samtools flagstat  ${sample}.rmdup.bam > ${sample}.rmdup.stat
	samtools view  -h  -f 2 -q 30    ${sample}.rmdup.bam   |grep -v chrM |samtools sort  -O bam  -@ 5 -o - > ${sample}.last.bam
	samtools index   ${sample}.last.bam 
	samtools flagstat  ${sample}.last.bam > ${sample}.last.stat 
	bedtools bamtobed -i ${sample}.last.bam  > ${sample}.bed
done
