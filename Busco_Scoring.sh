


##qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco_BP -q threaded.q -pe thread 8 -b y -cwd /home/jmattick/SNP_Scripts/Busco_Scoring.sh
##Make protein fasta of augustus output
if [[ 1 -eq 2 ]]
then
	/usr/local/packages/python-3.5/bin/python3 /usr/local/packages/augustus-3.3.1/scripts/getAnnoFastaFromJoingenes.py -s FILTER -g /autofs/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BrakerTemp/Braker/genome.fa -f /autofs/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BrakerTemp/Braker/augustus.hints.gtf -o /autofs/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BrakerTemp/Braker/augustus.hints.fa
fi

if [[ 1 -eq 2 ]]
then
	#/usr/local/packages/python-3.5/bin/python3 /usr/local/packages/augustus-3.3.1/scripts/getAnnoFastaFromJoingenes.py -s FILTER -g /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta -f /local/scratch/jmattick/LongReadMapping/DirectRNA.sorted.stringtie.gtf -o /local/scratch/jmattick/LongReadMapping/DirectRNA.sorted.stringtie.gtf
	#/usr/local/packages/cufflinks/bin/gffread /local/scratch/jmattick/LongReadMapping/DirectRNA.sorted.stringtie.gtf -g /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta
	#/usr/local/packages/bedtools2/bin/bedtools getfasta -fi /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta -bed /local/scratch/jmattick/LongReadMapping/DirectRNA.sorted.stringtie.gtf
#
	for i in $(cat /local/scratch/jmattick/LongReadMapping/DirectRNA.sorted.stringtie.gtf | grep -v "#" | awk -F '\t' '{print $9}' | sed 's/";.*//g' | sed 's/gene_id "//g' | sort | uniq)
	do
		echo $i
		qsub -P jhotopp-gcid-proj4b-filariasis -o /dev/null -e /dev/null -l mem_free=200M -b y -cwd /home/jmattick/SNP_Scripts/Busco_Scoring_mod.sh -i $i
	done
fi

if [[ 1 -eq 2 ]]
then
	cat /local/scratch/jmattick/LongReadMapping/*.fa > /local/scratch/jmattick/Stringtie.nuc.fa
	java -ea -Xmx10000m -cp /home/jmattick/SNP_Scripts/bbmap/current/ jgi.TranslateSixFrames in=/local/scratch/jmattick/Stringtie.nuc.fa out=/local/scratch/jmattick/Stringtie.AA.fa frames=1 overwrite=t

fi

if [[ 1 -eq 2 ]]
then
	rm /local/scratch/jmattick/Busco_BP.o*
	rm /local/scratch/jmattick/Busco_BP.e*
	rm -rf /local/scratch/jmattick/BUSCO_Scoring/
	mkdir /local/scratch/jmattick/BUSCO_Scoring/
	cp /autofs/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BrakerTemp/Braker/augustus.hints.fa.aa /local/scratch/jmattick/BUSCO_Scoring/Genes.aa.fa
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/python-3.5/lib/
	cp ~/config.ini /local/scratch/jmattick/BUSCO_Scoring/config.ini
	export BUSCO_CONFIG_FILE="/local/scratch/jmattick/BUSCO_Scoring/config.ini"
	cd /local/scratch/jmattick/BUSCO_Scoring/
	/usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/scripts/run_BUSCO.py -i /local/scratch/jmattick/BUSCO_Scoring/Genes.aa.fa -f -c 8 -l /home/jmattick/eukaryota_odb9/ -o BPahangi.Busco -m prot

	
fi

if [[ 1 -eq 2 ]]
then
	rm /local/scratch/jmattick/Busco_BP.o*
	rm /local/scratch/jmattick/Busco_BP.e*
	cp /local/aberdeen2rw/julie/JM_dir/GenomeSizes/Brugia_Malayi/brugia_malayi.PRJNA10729.WBPS11.protein.fa /local/scratch/jmattick/BUSCO_Scoring/BM.Genes.aa.fa
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/python-3.5/lib/
	cp ~/config.ini /local/scratch/jmattick/BUSCO_Scoring/config.ini
	export BUSCO_CONFIG_FILE="/local/scratch/jmattick/BUSCO_Scoring/config.ini"
	cd /local/scratch/jmattick/BUSCO_Scoring/
	/usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/scripts/run_BUSCO.py -i /local/scratch/jmattick/BUSCO_Scoring/BM.Genes.aa.fa -f -c 8 -l /home/jmattick/eukaryota_odb9/ -o BMalayi.Busco -m prot

fi


if [[ 1 -eq 2 ]]
then
	rm /local/scratch/jmattick/Busco_BP.o*
	rm /local/scratch/jmattick/Busco_BP.e*
	cp /local/aberdeen2rw/julie/JM_dir/GenomeSizes/Brugia_Pahangi/brugia_pahangi.PRJEB497.WBPS14.protein.fa /local/scratch/jmattick/BUSCO_Scoring/BP.old.Genes.aa.fa
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/python-3.5/lib/
	cp ~/config.ini /local/scratch/jmattick/BUSCO_Scoring/config.ini
	export BUSCO_CONFIG_FILE="/local/scratch/jmattick/BUSCO_Scoring/config.ini"
	cd /local/scratch/jmattick/BUSCO_Scoring/
	/usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/scripts/run_BUSCO.py -i /local/scratch/jmattick/BUSCO_Scoring/BP.old.Genes.aa.fa -f -c 8 -l /home/jmattick/eukaryota_odb9/ -o BPahangi.Old.Busco -m prot

fi

if [[ 1 -eq 2 ]]
then
	rm /local/scratch/jmattick/Busco_BP.o*
	rm /local/scratch/jmattick/Busco_BP.e*
	cp /local/scratch/jmattick/Stringtie.AA.fa /local/scratch/jmattick/BUSCO_Scoring/BP.stringtie.Genes.aa.fa
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/python-3.5/lib/
	cp ~/config.ini /local/scratch/jmattick/BUSCO_Scoring/config.ini
	export BUSCO_CONFIG_FILE="/local/scratch/jmattick/BUSCO_Scoring/config.ini"
	cd /local/scratch/jmattick/BUSCO_Scoring/
	/usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/scripts/run_BUSCO.py -i /local/scratch/jmattick/BUSCO_Scoring/BP.stringtie.Genes.aa.fa -f -c 8 -l /home/jmattick/eukaryota_odb9/ -o BPahangi.stringtie.Busco -m prot

fi

if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BUSCO_Finding/
	mkdir /local/scratch/jmattick/BUSCO_Finding/
	cd /local/scratch/jmattick/BUSCO_Finding/

	samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BM.Genes.aa.fa
	samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BP.old.Genes.aa.fa
	makeblastdb -in /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta -dbtype nucl

	cp ~/BMGenes_Busco.tsv /local/scratch/jmattick/BUSCO_Finding/
	cp ~/BP_oldGenes_Busco.tsv /local/scratch/jmattick/BUSCO_Finding/

	for i in $(cat /local/scratch/jmattick/BUSCO_Finding/BMGenes_Busco.tsv | awk '{print $2}' | sort | uniq)
	do
		mkdir /local/scratch/jmattick/BUSCO_Finding/$i/
		cd /local/scratch/jmattick/BUSCO_Finding/$i/
		BMGeneList=$(cat /local/scratch/jmattick/BUSCO_Finding/BMGenes_Busco.tsv | grep $i | awk '{print $1}' | tr '\n' ' ')
		samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BM.Genes.aa.fa $BMGeneList > /local/scratch/jmattick/BUSCO_Finding/$i/BMProteins.fa
		BPoldGeneList=$(cat /local/scratch/jmattick/BUSCO_Finding/BP_oldGenes_Busco.tsv | grep $i | awk '{print $1}' | tr '\n' ' ')
		samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BP.old.Genes.aa.fa $BPoldGeneList > /local/scratch/jmattick/BUSCO_Finding/$i/BPoldProteins.fa
		cat /local/scratch/jmattick/BUSCO_Finding/$i/BMProteins.fa /local/scratch/jmattick/BUSCO_Finding/$i/BPoldProteins.fa > /local/scratch/jmattick/BUSCO_Finding/$i/Proteins.fa
		qsub -P jhotopp-gcid-proj4b-filariasis -o /local/scratch/jmattick/BUSCO_Finding/$i/BPMatches.$i.blast -l mem_free=4G -q threaded.q -pe thread 12 -b y -cwd tblastn -query /local/scratch/jmattick/BUSCO_Finding/$i/Proteins.fa -db /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta -num_threads 12 -outfmt 6
	done

fi


if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BUSCO_Finding_Promer/
	mkdir /local/scratch/jmattick/BUSCO_Finding_Promer/
	cd /local/scratch/jmattick/BUSCO_Finding_Promer/

	samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BM.Genes.nuc.fa
	samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BP.old.Genes.nuc.fa

	cp ~/BMGenes_Busco.tsv /local/scratch/jmattick/BUSCO_Finding_Promer/
	cp ~/BP_oldGenes_Busco.tsv /local/scratch/jmattick/BUSCO_Finding_Promer/

	for i in $(cat /local/scratch/jmattick/BUSCO_Finding_Promer/BMGenes_Busco.tsv | awk '{print $2}' | sort | uniq)
	do
		mkdir /local/scratch/jmattick/BUSCO_Finding_Promer/$i/
		cd /local/scratch/jmattick/BUSCO_Finding_Promer/$i/
		BMGeneList=$(cat /local/scratch/jmattick/BUSCO_Finding_Promer/BMGenes_Busco.tsv | grep $i | awk '{print $1}' | tr '\n' ' ')
		samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BM.Genes.nuc.fa $BMGeneList > /local/scratch/jmattick/BUSCO_Finding_Promer/$i/BMProteins.fa
		echo "Empty?"
		head -2 /local/scratch/jmattick/BUSCO_Finding_Promer/$i/BMProteins.fa
		BPoldGeneList=$(cat /local/scratch/jmattick/BUSCO_Finding_Promer/BP_oldGenes_Busco.tsv | grep $i | awk '{print $1}' | tr '\n' ' ')
		samtools faidx /local/scratch/jmattick/BUSCO_Scoring/BP.old.Genes.nuc.fa $BPoldGeneList > /local/scratch/jmattick/BUSCO_Finding_Promer/$i/BPoldProteins.fa
		echo "Empty?"
		head -2 /local/scratch/jmattick/BUSCO_Finding_Promer/$i/BPoldProteins.fa
		cat /local/scratch/jmattick/BUSCO_Finding_Promer/$i/BMProteins.fa /local/scratch/jmattick/BUSCO_Finding_Promer/$i/BPoldProteins.fa > /local/scratch/jmattick/BUSCO_Finding_Promer/$i/Proteins.fa
		qsub -P jhotopp-gcid-proj4b-filariasis  -l mem_free=10G -b y -cwd /usr/local/packages/mummer/promer -p FindingGenes.$i /local/scratch/jmattick/BUSCO_Finding_Promer/$i/Proteins.fa /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta
	done

fi


if [[ 1 -eq 2 ]]
then
	for i in $(cat /local/scratch/jmattick/BUSCO_Finding_Promer/BMGenes_Busco.tsv | awk '{print $2}' | sort | uniq)
	do
			/usr/local/packages/mummer/show-coords -qlTHb /local/scratch/jmattick/BUSCO_Finding_Promer/$i/FindingGenes.$i.delta > /local/scratch/jmattick/BUSCO_Finding_Promer/$i/FindingGenes.$i.coords
	done
fi
