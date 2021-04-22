ls ~/login123/ATAC-seq/bowite2_result/last_bam/*.last.bam|while read id
do 
    nohup /ifs1/User/shangyao/miniconda3/bin/bamCoverage -p 15 --binSize 10 --normalizeUsing CPM -b $id -o ${id%%.*}.last.bw >${id%%.*}.log 2>&1 &
done
