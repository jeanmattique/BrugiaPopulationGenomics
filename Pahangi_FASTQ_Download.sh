rm -rf /local/scratch/jmattick/PahangiFASTQ/
mkdir /local/scratch/jmattick/PahangiFASTQ/

for i in $(cat /home/jmattick/SNP_Scripts/BPahangiSRA.tsv | awk '{print $1}' | sort | uniq)
do
	mkdir /local/scratch/jmattick/PahangiFASTQ/$i
done



while read line
do
	dir=$(echo $line | awk '{print $1}')
	SRR=$(echo $line | awk '{print $2}')
	qsub -P jhotopp-gcid-proj4b-filariasis -o /local/scratch/jmattick/PahangiFASTQ/$SRR.out -e /local/scratch/jmattick/PahangiFASTQ/$SRR.error -l mem_free=2G -cwd /home/jmattick/SNP_Scripts/Pahangi_FASTQ_Download_mod.sh -i $dir -j $SRR
	#cd /local/scratch/jmattick/PahangiFASTQ/$dir
	#/usr/local/packages/sratoolkit-2.9.0/bin/fastq-dump --split-files $SRR
done < /home/jmattick/SNP_Scripts/BPahangiSRA.tsv