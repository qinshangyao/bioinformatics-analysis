for i in `ls *.gz`;
do 
    echo $i
    fastqc $i -o ../fastqc_result -q
done
