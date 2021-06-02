for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-gene" ]] 
		then
		var=$(($var+1))
		gene=${!var}
	elif [[ "${!var}" = "-chr" ]] 
		then
		var=$(($var+1))
		chr=${!var}
	elif [[ "${!var}" = "-same" ]] 
		then
		var=$(($var+1))
		same=${!var}
	else
		exit 1
	fi
done

BPPath=$(echo "/local/scratch/jmattick/BUSCO_Predictions/BPahangi/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".faa")
BMPath=$(echo "/local/scratch/jmattick/BUSCO_Predictions/BMalayi/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".faa")
WBPath=$(echo "/local/scratch/jmattick/BUSCO_Predictions/WBancrofti/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".faa")
BTPath=$(echo "/local/scratch/jmattick/BUSCO_Predictions/BTimori/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".faa")
OVPath=$(echo "/local/scratch/jmattick/BUSCO_Predictions/OVolvulus/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".faa")

cat $BPPath $BMPath $WBPath $BTPath $OVPath | awk '{print $1}' > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.faa

/usr/local/packages/clustalo/bin/clustalo -i /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.faa --distmat-out /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.tsv -o /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.aligned.faa --percent-id --force --full

maxid=$(cat /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.tsv | tail -n+2 | awk '{print $2"\n"$3"\n"$4"\n"$5"\n"$6}' | awk 'NR != 1 && NR != 7 && NR != 13 && NR != 19 && NR != 25 {print $0}' | sort -k1,1nr | head -1)
minid=$(cat /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.tsv | tail -n+2 | awk '{print $2"\n"$3"\n"$4"\n"$5"\n"$6}' | awk 'NR != 1 && NR != 7 && NR != 13 && NR != 19 && NR != 25 {print $0}' | sort -k1,1n | head -1)

minlen=$(/home/jdhotopp/bin/residues.pl /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.faa | awk '{print $2}' | sort -k1,1n | head -1)
ratio=$(/home/jdhotopp/bin/residues.pl /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.faa | awk '{print $2}' | sort -k1,1nr | head -1 | awk -v minlen="$minlen" '{print 100*(minlen/$1)}')
echo $ratio

echo -e $gene"\t"$chr"\t"$same"\t"$maxid"\t"$minid"\t"$ratio > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.output



BPPathNuc=$(echo "/local/scratch/jmattick/BUSCO_Predictions/BPahangi/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".fna")
BMPathNuc=$(echo "/local/scratch/jmattick/BUSCO_Predictions/BMalayi/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".fna")
WBPathNuc=$(echo "/local/scratch/jmattick/BUSCO_Predictions/WBancrofti/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".fna")
BTPathNuc=$(echo "/local/scratch/jmattick/BUSCO_Predictions/BTimori/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".fna")
OVPathNuc=$(echo "/local/scratch/jmattick/BUSCO_Predictions/OVolvulus/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/"$gene".fna")

cat $BPPathNuc $BMPathNuc $WBPathNuc $BTPathNuc $OVPathNuc | awk '{print $1}' > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.fna

perl /home/jmattick/SNP_Scripts/TranslatorX.pl -i /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.fna -o /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned

BPContig=$(cat $BPPathNuc | awk '{print $1}' | grep ">" | sed 's/>//g')
BMContig=$(cat $BMPathNuc | awk '{print $1}' | grep ">" | sed 's/>//g')
WBContig=$(cat $WBPathNuc | awk '{print $1}' | grep ">" | sed 's/>//g')
BTContig=$(cat $BTPathNuc | awk '{print $1}' | grep ">" | sed 's/>//g')
OVContig=$(cat $OVPathNuc | awk '{print $1}' | grep ">" | sed 's/>//g')

samtools faidx /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.fasta

samtools faidx /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.fasta $BPContig > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.BP.fasta
samtools faidx /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.fasta $BMContig > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.BM.fasta
samtools faidx /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.fasta $WBContig > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.WB.fasta
samtools faidx /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.fasta $BTContig > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.BT.fasta
samtools faidx /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.fasta $OVContig > /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/$gene.AAaligned.nt_ali.OV.fasta



