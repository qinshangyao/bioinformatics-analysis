bed=~/ep_opc/olig2_binding_site.bed
for id in *.bw;
do
echo $id
file=$(basename $id)
sample=${file%%.*}
echo $sample

/ifs1/User/shangyao/miniconda3/bin/computeMatrix reference-point --referencePoint TSS -p 30 -b 2000 -a 2000 -R $bed -S $id --skipZeros -o matrix_tss2k_${sample}.gz --outFileSortedRegions region_tss2k_${sample}.bed

plotHeatmap -m matrix_tss2k_${sample}.gz -out ${sample}_heatmap_2k.png
plotProfile -m matrix_tss2k_${sample}.gz -out ${sample}_profile_2k.png

done
