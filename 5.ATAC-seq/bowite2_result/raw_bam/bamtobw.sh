ls *.bam | while read id
do
    nohup  ~/downlaod/dt/deepTools/bin/bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.rm.bw &
done
