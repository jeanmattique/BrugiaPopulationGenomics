for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-i" ]] 
		then
		var=$(($var+1))
		i=${!var}
		echo "Variable i is "$i
	else
		echo "Unrecognized input, exitting"
		exit 1
	fi
done


echo ">"$i > /local/scratch/jmattick/LongReadMapping/Temp.$i.fa
cat /local/scratch/jmattick/LongReadMapping/DirectRNA.sorted.stringtie.gtf | grep -v "#" | grep -w $i | awk '$3 == "exon" {print $0}' > /local/scratch/jmattick/LongReadMapping/Temp.$i.gtf
direction=$(cat /local/scratch/jmattick/LongReadMapping/Temp.$i.gtf | awk '{print $7}' | sort | uniq)
if [ $direction == "+" ]
then
	/usr/local/packages/bedtools2/bin/bedtools getfasta -fi /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta -bed /local/scratch/jmattick/LongReadMapping/Temp.$i.gtf | grep -v ">" | tr -d '\n' >> /local/scratch/jmattick/LongReadMapping/Temp.$i.fa
else
	/usr/local/packages/bedtools2/bin/bedtools getfasta -fi /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta -bed /local/scratch/jmattick/LongReadMapping/Temp.$i.gtf | grep -v ">" | tr -d '\n' | rev | tr ACGTacgt TGCAtgca >> /local/scratch/jmattick/LongReadMapping/Temp.$i.fa
fi
rm /local/scratch/jmattick/LongReadMapping/Temp.$i.gtf
samtools faidx /local/scratch/jmattick/LongReadMapping/Temp.$i.fa
samtools faidx /local/scratch/jmattick/LongReadMapping/Temp.$i.fa $i > /local/scratch/jmattick/LongReadMapping/$i.fa
rm /local/scratch/jmattick/LongReadMapping/Temp.$i.fa
rm /local/scratch/jmattick/LongReadMapping/Temp.$i.fa.fai
