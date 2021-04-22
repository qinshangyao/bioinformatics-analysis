for i in $(seq 12 23);
do
    samtools view -S SRR115109${i}.sam -b > SRR115109${i}.bam
    samtools sort SRR115109${i}.bam -o SRR115109${i}_sorted.bam
    samtools index SRR115109${i}_sorted.bam
done
