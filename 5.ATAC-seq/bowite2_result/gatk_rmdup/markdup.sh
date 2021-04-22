for i in `ls *.raw.bam | cut -d "." -f 1`
do
    echo "processing $i"
    nohup gatk MarkDuplicates -I ${i}.raw.bam -O ${i}.picard.rmdup.bam --REMOVE_SEQUENCING_DUPLICATES true -M ${i}.picard.log >${i}.nohup.log 2>&1 &
done
