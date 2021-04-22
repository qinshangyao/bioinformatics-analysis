for fn in `cat sampleNames_12`
do 
    echo "Processing sample ${fn}"
    salmon quant -i Mus_musculus_ensembl_index/ -l A -1 ${fn}_1.fastq.gz -2 ${fn}_2.fastq.gz -p 8 --validateMappings -o quants/${fn}_quant
done 
