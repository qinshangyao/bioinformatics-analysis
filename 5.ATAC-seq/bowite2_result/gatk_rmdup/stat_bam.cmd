for i in `ls *.picard.rmdup.bam`
do
    samtools flagstat $i > ${i}.stat
done
