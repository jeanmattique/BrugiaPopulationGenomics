if [[ 1 -eq 2 ]]
then
	#rm -rf /local/scratch/jmattick/BUSCO_Predictions/
	#mkdir /local/scratch/jmattick/BUSCO_Predictions/

	cd /local/scratch/jmattick/BUSCO_Predictions/
	#export BUSCO_CONFIG_FILE="/usr/local/packages/busco/config/config.ini"
	cat /usr/local/packages/busco/config/config.ini | sed 's/\/ncbi-blast-2.2.31+/\/usr\/local\/packages\/ncbi-blast+/g' | sed 's/\/augustus/\/usr\/local\/packages\/augustus/g' > /local/scratch/jmattick/BUSCO_Predictions/config.ini
	export AUGUSTUS_CONFIG_PATH=/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/AUGUSTUS_Config/augustus-3.3.3/config/
	export AUGUSTUS_BIN_PATH=/usr/local/packages/augustus-3.3.3/bin/
	export AUGUSTUS_SCRIPTS_PATH=/usr/local/packages/augustus-3.3.3/scripts/

	export BUSCO_CONFIG_FILE="/local/scratch/jmattick/BUSCO_Predictions/config.ini"
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/python-3.5/lib/

	fastapath=$(echo "/local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.fasta")
	#qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o BPahangi -l nematoda_odb10 --update-data > /local/scratch/jmattick/BUSCO_Predictions/nematoda.update.txt
	#list=$(cat /local/scratch/jmattick/BUSCO_Predictions/nematoda.update.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o BPahangi -l nematoda_odb10


	fastapath=$(echo "/local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa")

	#qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o BMalayi -l nematoda_odb10
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o BMalayi -l nematoda_odb10

	fastapath=$(echo "/local/scratch/jmattick/wuchereria_bancrofti.PRJEB536.WBPS15.genomic.fa")

	#qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o WBancrofti -l nematoda_odb10
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o WBancrofti -l nematoda_odb10


	fastapath=$(echo "/local/scratch/jmattick/brugia_timori.PRJEB4663.WBPS15.genomic.fa")

	#qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o BTimori -l nematoda_odb10
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o BTimori -l nematoda_odb10

	fastapath=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/ovolvulus.genome.fa")

	#qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o OVolvulus -l nematoda_odb10
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -N Busco -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/python-3.5/bin/python3.5 /usr/local/packages/busco/bin/busco -m genome -i $fastapath -f -o OVolvulus -l nematoda_odb10


fi

if [[ 1 -eq 2 ]]
then
	rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/


	cp /local/scratch/jmattick/BUSCO_Predictions/BPahangi/run_nematoda_odb10/full_table.tsv /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/BPahangi.fulltable.tsv
	cp /local/scratch/jmattick/BUSCO_Predictions/BMalayi/run_nematoda_odb10/full_table.tsv /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/BMalayi.fulltable.tsv
	cp /local/scratch/jmattick/BUSCO_Predictions/WBancrofti/run_nematoda_odb10/full_table.tsv /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/WBancrofti.fulltable.tsv
	cp /local/scratch/jmattick/BUSCO_Predictions/BTimori/run_nematoda_odb10/full_table.tsv /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/BTimori.fulltable.tsv
	cp /local/scratch/jmattick/BUSCO_Predictions/OVolvulus/run_nematoda_odb10/full_table.tsv /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/OVolvulus.fulltable.tsv



	##R: Final_PopGenome_Figures.R
	##Busco Multispecies section
	##Output: /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/ROutput.busco.multispecies.tsv
fi

if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/
	mkdir /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/

	cd /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/

	while read line
	do
		genename=$(echo $line | awk '{print $1}')
		echo $genename
		chr=$(echo $line | awk '{print $7}')
		same=$(echo $line | awk '{print $9}')
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=1G -b y -cwd /home/jmattick/SNP_Scripts/BuscoGeneEval.sh -gene $genename -chr $chr -same $same
	done < /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/ROutput.busco.multispecies.tsv
fi

if [[ 1 -eq 2 ]]
then
	cat /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/*.output > /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv
	rm -rf /local/scratch/jmattick/BUSCO_Predictions/IndAlignments/
	mkdir /local/scratch/jmattick/BUSCO_Predictions/IndAlignments/
	cp /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/*.AAaligned.nt_ali.B*.fasta /local/scratch/jmattick/BUSCO_Predictions/IndAlignments/
	cp /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/*.AAaligned.nt_ali.WB.fasta /local/scratch/jmattick/BUSCO_Predictions/IndAlignments/
	cp /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/*.AAaligned.nt_ali.OV.fasta /local/scratch/jmattick/BUSCO_Predictions/IndAlignments/

	rm -rf /local/scratch/jmattick/BUSCO_Predictions/GeneAssessment/
fi


###Filter gene list for 85% minimum match across all 4 species and all genes are within 90% of length
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/
	mkdir /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/

	cd /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/

	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.BP.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BPahangi.ChrX.genelist
	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep -v "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.BP.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BPahangi.autosome.genelist

	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.BM.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BMalayi.ChrX.genelist
	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep -v "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.BM.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BMalayi.autosome.genelist

	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.WB.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/WBancrofti.ChrX.genelist
	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep -v "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.WB.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/WBancrofti.autosome.genelist

	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.BT.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BTimori.ChrX.genelist
	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep -v "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.BT.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BTimori.autosome.genelist


	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.OV.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/OVolvulus.ChrX.genelist
	cat /local/scratch/jmattick/BUSCO_Predictions/AlignmentQual.tsv | grep "Same" | grep -v "ChrX" | awk '$5 > 85 && $6 > 90 {print "/local/scratch/jmattick/BUSCO_Predictions/IndAlignments/"$1".AAaligned.nt_ali.OV.fasta"}' > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/OVolvulus.autosome.genelist

	echo ">BPahangi" > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BPahangi.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	echo -e "\n>BMalayi" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BMalayi.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	echo -e "\n>WBancrofti" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/WBancrofti.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	echo -e "\n>BTimori" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BTimori.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	echo -e "\n>OVolvulus" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/OVolvulus.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna


	echo ">BPahangi" > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BPahangi.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	echo -e "\n>BMalayi" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BMalayi.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	echo -e "\n>WBancrofti" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/WBancrofti.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	echo -e "\n>BTimori" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BTimori.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	echo -e "\n>OVolvulus" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/OVolvulus.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna


	/home/jmattick/SNP_Scripts/convertFasta2Phylip.sh /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.phy
	/home/jmattick/SNP_Scripts/convertFasta2Phylip.sh /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.phy

	cd /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/
	#/usr/local/packages/raxml/bin/raxmlHPC -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMA -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.phy -n RaxML.XChrom
	#/usr/local/packages/raxml/bin/raxmlHPC -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMA -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.phy -n RaxML.autosome
	/usr/local/packages/iqtree/bin/iqtree -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.phy -nt 4 -bb 10000 -redo
	/usr/local/packages/iqtree/bin/iqtree -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.phy -nt 4 -bb 10000 -redo
	##IQ Tree? Model Testing. ProtTest

	cp /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/RAxML_bipartitionsBranchLabels.RaxML.XChrom /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/RAxML_bipartitionsBranchLabels.RaxML.XChrom
	cp /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/RAxML_bipartitionsBranchLabels.RaxML.autosome /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/RAxML_bipartitionsBranchLabels.RaxML.autosome

	##Wol phylogeny -- ?
	##Mito phylogeny -- core ungapped alignment

	##Freq of w bancrofti mpileup alt bases
fi

###DO ABOVE EXCEPT WITHOUT BRUGIA TIMORI
if [[ 1 -eq 2 ]]
then
	cd /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/

	echo ">BPahangi" > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BPahangi.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna
	echo -e "\n>BMalayi" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BMalayi.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna
	echo -e "\n>WBancrofti" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/WBancrofti.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna
	echo -e "\n>OVolvulus" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/OVolvulus.ChrX.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.noBT.fna


	echo ">BPahangi" > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BPahangi.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna
	echo -e "\n>BMalayi" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/BMalayi.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna
	echo -e "\n>WBancrofti" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/WBancrofti.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna
	echo -e "\n>OVolvulus" >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna
	grep -v '^#' /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/OVolvulus.autosome.genelist | xargs cat | grep -v ">" | tr -d '\n' | fold -w 80 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.noBT.fna

fi

###Identify and assemble B. timori mitochondria and wolbachia
if [[ 1 -eq 2 ]]
then

	######WOLBACHIA
	outdir=$(echo "/local/scratch/jmattick/BUSCO_Predictions/Timori_Wolbachia/")
	rm -rf $outdir
	mkdir $outdir
	#mv ~/brugia_timori.PRJEB4663.WBPS14.genomic.fa.gz /local/scratch/jmattick/
	#gzip -d /local/scratch/jmattick/brugia_timori.PRJEB4663.WBPS14.genomic.fa.gz
	fastapath=$(echo "/local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.wol.fa")
	/usr/local/packages/bwa/bin/bwa index $fastapath
	/usr/local/packages/samtools/bin/samtools faidx $fastapath
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath


	fastq1=$(echo "/local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_1.fastq")
	fastq2=$(echo "/local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_2.fastq")

	j="wBT"
	qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -k 23 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

	list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
		
	list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPDedup.txt

	list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools  depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam



	#####MITOCHONDRIA
	outdir=$(echo "/local/scratch/jmattick/BUSCO_Predictions/Timori_Mitochondria/")
	rm -rf $outdir
	mkdir $outdir
	#mv ~/brugia_timori.PRJEB4663.WBPS14.genomic.fa.gz /local/scratch/jmattick/
	#gzip -d /local/scratch/jmattick/brugia_timori.PRJEB4663.WBPS14.genomic.fa.gz
	fastapath=$(echo "/local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.mt.fa")
	/usr/local/packages/bwa/bin/bwa index $fastapath
	/usr/local/packages/samtools/bin/samtools faidx $fastapath
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath


	fastq1=$(echo "/local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_1.fastq")
	fastq2=$(echo "/local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_2.fastq")

	j="mBT"
	qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -k 23 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

	list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
		
	list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPDedup.txt

	list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools  depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

fi

if [[ 1 -eq 2 ]]
then

	######WOLBACHIA
	outdir=$(echo "/local/scratch/jmattick/BUSCO_Predictions/Timori_Assembly/")
	rm -rf $outdir
	mkdir $outdir
	samtools view -F4 /local/scratch/jmattick/BUSCO_Predictions/Timori_Wolbachia/wBT.sorted.bam | awk '{print $1}' > $outdir/wBt.reads
	samtools view -F4 /local/scratch/jmattick/BUSCO_Predictions/Timori_Mitochondria/mBT.sorted.bam | awk '{print $1}' > $outdir/mBt.reads
	cat $outdir/wBt.reads | sort | uniq > $outdir/wBt.uniq.reads
	cat $outdir/mBt.reads | sort | uniq > $outdir/mBt.uniq.reads

	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar FilterSamReads I=/local/scratch/jmattick/BUSCO_Predictions/Timori_Wolbachia/wBT.sorted.bam O=$outdir/wBt.filtered.bam FILTER=includeReadList READ_LIST_FILE=$outdir/wBt.uniq.reads
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar FilterSamReads I=/local/scratch/jmattick/BUSCO_Predictions/Timori_Mitochondria/mBT.sorted.bam O=$outdir/mBT.filtered.bam FILTER=includeReadList READ_LIST_FILE=$outdir/mBt.uniq.reads
	cd $outdir
	/usr/local/packages/samtools/bin/samtools bam2fq -1 wBt.filtered_1.fastq -2 wBt.filtered_2.fastq wBt.filtered.bam
	/usr/local/packages/samtools/bin/samtools bam2fq -1 mBT.filtered_1.fastq -2 mBT.filtered_2.fastq mBT.filtered.bam

	PATH=/home/jmattick/.local/bin:$PATH
	rm -rf $outdir/GetOrganelleAll/
	mkdir $outdir/GetOrganelleAll/

	qsub -V -P jhotopp-gcid-proj4b-filariasis -N GetOrganelle -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd get_organelle_from_reads.py -1 /local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_1.fastq -2 /local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_2.fastq --reduce-reads-for-coverage inf -w 0.7 -o $outdir/GetOrganelleAll/ -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.mt.fa -t 32 -F animal_mt

	rm -rf $outdir/GetOrganelleFiltered/
	mkdir $outdir/GetOrganelleFiltered/

	qsub -V -P jhotopp-gcid-proj4b-filariasis -N GetOrganelle -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd get_organelle_from_reads.py -1 $outdir/mBT.filtered_1.fastq -2 $outdir/mBT.filtered_2.fastq --reduce-reads-for-coverage inf -w 0.7 -o $outdir/GetOrganelleFiltered/ -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.mt.fa -t 32 -F animal_mt



	use unicycler
	rm -rf $outdir/Unicycler/
	mkdir $outdir/Unicycler/
	qsub -P jhotopp-gcid-proj4b-filariasis -V -l mem_free=8G -q threaded.q -pe thread 8 -b y -cwd unicycler --mode bold --largest_component -1 $outdir/wBt.filtered_1.fastq -2 $outdir/wBt.filtered_2.fastq -o $outdir/Unicycler/ -t 8 --pilon_path /usr/local/packages/pilon-1.22/pilon-1.22.jar

	rm -rf $outdir/mt_Unicycler/
	mkdir $outdir/mt_Unicycler/
	qsub -P jhotopp-gcid-proj4b-filariasis -V -l mem_free=8G -q threaded.q -pe thread 8 -b y -cwd unicycler --mode bold --largest_component -1 $outdir/mBT.filtered_1.fastq -2 $outdir/mBT.filtered_2.fastq -o $outdir/mt_Unicycler/ -t 8 --pilon_path /usr/local/packages/pilon-1.22/pilon-1.22.jar

	cat $outdir/GetOrganelleFiltered/animal_mt.K85.contigs.graph1.1.path_sequence.fasta | sed 's/--.*//g' > BT.filteredassembly.fa

	nucmer -p BTMitoMap /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.mt.fa $outdir/BT.filteredassembly.fa
	/usr/local/packages/mummer/show-coords -qlTHb BTMitoMap.delta > BTMitoMap.coords

	samtools faidx $outdir/BT.filteredassembly.fa contig_1 > $outdir/BT.filteredassembly.contig1.fa
	samtools faidx $outdir/BT.filteredassembly.fa contig_2 > $outdir/BT.filteredassembly.contig2.fa

	nucmer -p BTMitoSelf $outdir/BT.filteredassembly.contig1.fa $outdir/BT.filteredassembly.contig2.fa
	/usr/local/packages/mummer/show-coords -qlTHb BTMitoSelf.delta > BTMitoSelf.coords


	###Rotation Time
	samtools faidx $outdir/BT.filteredassembly.fa contig_2:2962-7275 > $outdir/BT.rotated.trimmed.fa
	samtools faidx $outdir/BT.filteredassembly.fa contig_1 >> $outdir/BT.rotated.trimmed.fa
	samtools faidx $outdir/BT.filteredassembly.fa contig_2:1-2961 >> $outdir/BT.rotated.trimmed.fa

	echo ">BT_Mitochondria" > $outdir/BT.rotated.trimmed.temp.fa
	cat $outdir/BT.rotated.trimmed.fa | grep -v ">" | tr -d '\n' >> $outdir/BT.rotated.trimmed.temp.fa
	samtools faidx $outdir/BT.rotated.trimmed.temp.fa BT_Mitochondria > $outdir/BT.rotated.trimmed.final.fa

	nucmer -p BTvsBMFinal /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.mt.fa $outdir/BT.rotated.trimmed.final.fa
	/usr/local/packages/mummer/show-coords -qlTH BTvsBMFinal.delta > BTvsBMFinal.coords

fi


if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/
	mkdir /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/

	mtcontig=$(echo "gi|23307675|gb|AF538716.1|")
	samtools faidx /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa $mtcontig > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.mt.fa
	wolcontig=$(echo "Bm_006")

	samtools faidx /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa $wolcontig > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.wol.fa

	cd /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/

	nucmer -p BTMitoSearch /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.mt.fa /local/scratch/jmattick/brugia_timori.PRJEB4663.WBPS14.genomic.fa
	nucmer -p OVWolSearch /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BM.wol.fa /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/ovolvulus.genome.fa

	/usr/local/packages/mummer/show-coords -qlTHb BTMitoSearch.delta > BTMitoSearch.coords
	/usr/local/packages/mummer/show-coords -qlTHb OVWolSearch.delta > OVWolSearch.coords

	echo ">wWb_scaffolded" > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/wWb.scaffold.fa
	cat GCF_002204235.2_ASM220423v2_genomic.fna | grep -v ">" | awk '{print $1}' | tr -d '\n' | fold -w 70 >> /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/wWb.scaffold.fa
	samtools faidx /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/wWb.scaffold.fa wWb_scaffolded > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/wWb.scaffold.fixed.fa

	mtcontig=$(echo "OVOC_MITOCHONDRIAL")
	samtools faidx /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/ovolvulus.genome.fa $mtcontig > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/OV.mt.fa

	mtcontig=$(echo "BP_Mito_unicyclerAssembly_pilon")
	samtools faidx /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.fasta $mtcontig > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BP.mt.fa
	wolcontig=$(echo "wBP_unicyclerAssembly_pilon")

	samtools faidx /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta $wolcontig > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/BP.wol.fa

	cat BM.wol.fa BP.wol.fa wWb.scaffold.fixed.fa wOv.fasta > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Wolbachia.multifasta.fa
	###BT Added
	cp /local/scratch/jmattick/BUSCO_Predictions/Timori_Assembly/BT.rotated.trimmed.final.fa /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/
	cat BM.mt.fa ~/BT.genbank.fa BP.mt.fa WB_mt.fasta OV.mt.fa > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Mitochondria.multifasta.fa
	#cat BM.mt.fa BP.mt.fa WB_mt.fasta OV.mt.fa > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Mitochondria.multifasta.fa

	rm -rf ~/Filarial_MT_FNA
	mkdir ~/Filarial_MT_FNA


	samtools faidx BM.mt.fa gi\|23307675\|gb\|AF538716.1\| | sed 's/gi.*/BM_mt/g' > ~/Filarial_MT_FNA/BM.mt.fa
	samtools faidx ~/BT.genbank.fa AP017686.1 | sed 's/AP017686.1/BT_mt/g' > ~/Filarial_MT_FNA/BT.mt.fa
	samtools faidx BP.mt.fa BP_Mito_unicyclerAssembly_pilon | sed 's/BP_Mito_unicyclerAssembly_pilon/BP_mt/g' > ~/Filarial_MT_FNA/BP.mt.fa
	samtools faidx WB_mt.fasta NC_016186.1 | sed 's/NC_016186.1/WB_mt/g' > ~/Filarial_MT_FNA/WB.mt.fa
	samtools faidx OV.mt.fa OVOC_MITOCHONDRIAL | sed 's/OVOC_MITOCHONDRIAL/OV_mt/g' > ~/Filarial_MT_FNA/OV.mt.fa
	cat ~/Filarial_MT_FNA/*.fa > /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Mitochondria.multifasta.fixed.fa

	FNA="/local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Wolbachia.multifasta.fa"
	/usr/local/packages/mafft/bin/mafft "$FNA" > "$(echo "$FNA" | sed "s/[.]fa*/.aligned.fa/g")"
	FNA="/local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Mitochondria.multifasta.fixed.fa"
	/usr/local/packages/mafft/bin/mafft "$FNA" > "$(echo "$FNA" | sed "s/Mitochondria.multifasta.fixed.fa/Mitochondria.multifasta.fixed.aligned.fa/g")"


	/usr/local/packages/iqtree/bin/iqtree -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Mitochondria.multifasta.fixed.aligned.fa -nt 4 -bb 1000 -redo
	cp /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA_SubGenomes/Mitochondria.multifasta.fixed.aligned.fa.treefile ~/
	/usr/local/packages/iqtree/bin/iqtree -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna -nt 4 -bb 1000 -redo
	cp /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.Xchrom.fna.treefile ~/

	/usr/local/packages/iqtree/bin/iqtree -s /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna -nt 4 -bb 1000 -redo
	cp /local/scratch/jmattick/BUSCO_Predictions/MultiFASTA/MultiAlignment.autosome.fna.treefile ~/

	mkdir ~/Filarial_MT_FNA

fi
