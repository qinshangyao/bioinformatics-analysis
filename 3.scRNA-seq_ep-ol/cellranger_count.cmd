ref=~/login123/scRNA-seq/olig2_inj/ref/Mus_cellranger_index/Mus_GRCm38_tdtomato
fastqs=/ifs1/User/shangyao/login123/scRNA-seq/olig2_inj/fastq

for id in `cat srrfile`
do
    nohup cellranger count --id=${id} --transcriptome=${ref} --fastqs=${fastqs} --sample=${id} --expect-cells=1000 --localcores=8 --localmem=32 >${id}.countlog 2>&1 &
done
