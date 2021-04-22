index=/ifs1/User/shangyao/ncbi/public/sra/danRer_genome/danRer_salmon_index
cd /ifs1/User/shangyao/ncbi/public/sra/clean_data
for fn in *.fq.gz
do
    prefix=$(ls $fn| cut -d "_" -f 1)
    echo ${prefix}
    ~/downlaod/salmon-latest_linux_x86_64/bin/salmon quant -i ${index} -l A -r $fn -p 30 --validateMappings -o transcripts_quant/${prefix}
done
