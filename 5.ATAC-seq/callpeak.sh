for i in `cat prefix`
do
    nohup macs2 callpeak -t $i.bed -n $i -g mm --nomodel --shift 75 --extsize 150 --outdir ../../callpeak_result/ > ${i}.callpeak.log 2>&1 &
done
