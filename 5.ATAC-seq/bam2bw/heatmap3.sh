bed=~/ep_opc/olig2_binding_site.bed
for id in *.bw;
do
echo $id
file=$(basename $id)
sample=${file%%.*}
echo $sample

/ifs1/User/shangyao/miniconda3/bin/computeMatrix scale-regions  -p 30 -m 500 -R $bed -S $id --skipZeros -o matrix_tss2k_${sample}.gz --outFileSortedRegions region_tss2k_${sample}.bed

plotHeatmap --colorList white,yellow,red -m matrix_tss2k_${sample}.gz -out ${sample}_heatmap_2k.png -min 0.2 -max 1.6

#plotProfile -m matrix_tss2k_${sample}.gz -out ${sample}_profile_2k.png

done
