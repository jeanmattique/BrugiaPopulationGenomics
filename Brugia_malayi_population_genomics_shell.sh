if [[ 1 -eq 2 ]]
then
	rm -rf /local/aberdeen2rw/julie/JM_dir/BMalayiClinical/
	mkdir /local/aberdeen2rw/julie/JM_dir/BMalayiClinical/
fi


###Clinical Sample Mapping
if [[ 1 -eq 2 ]]
then
	#fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa)
	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	#/usr/local/packages/bwa/bin/bwa index $fastapath
	#rm $dictpath
	#/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

	outdir=$(echo "/local/scratch/jmattick/BMalayiMalaysia_unmasked")
	#outdir=$(echo "/local/scratch/jmattick/BMalayiMalaysia")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for j in $( \ls /local/projects/EBMAL/ | grep "Malaysia");
	do
		fastq1=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz)
		fastq2=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R2.fastq.gz)

		qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

		list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
		
		list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

		list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools  depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

		#qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk AddOrReplaceReadGroups --INPUT $outdir/$j.sorted.dedup.bam --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM $j --OUTPUT $outdir/$j.sorted.dedup.readgroup.bam --CREATE_INDEX true

		##Works. Why does BWA fuck it up?
		##samtools view -H $outdir/$j.sorted.dedup.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' |  samtools reheader - $outdir/$j.sorted.dedup.bam > $outdir/$j.sorted.dedup.reheader.bam

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf
		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --sample-ploidy 1 --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.haploid.vcf

	done
fi



##Mapping Stats
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/BMalayiStats/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa)


	for j in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Bams/ | grep "bam");
	do

		Sample=$(echo $j | awk -F '.' '{print $1".flagstat"}')
		qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$Sample -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Bams/$j
	done
	cp /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Bams/*.metrics $outdir


fi

##Build Sequencing Stats Table
if [[ 1 -eq 2 ]]
then
	echo -e "Sample\tTotalReads\tMappedReads\tPercentMapped\tPercentDuplicates\tGenomeCoverage" > ~/BM.Mapping.Summary.tsv
	GenomeSize=90163045
	for j in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Bams/ | grep "bam");
	do
		Sample=$(echo $j | awk -F '.' '{print $1}')
		metrics=$(echo $j | sed 's/bam/metrics/g')
		totalreads=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiStats/$Sample.flagstat | head -1 | sed 's/ +.*//g')
		totalmapped=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiStats/$Sample.flagstat | grep "mapped (" | sed 's/ +.*//g')
		percmapped=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiStats/$Sample.flagstat | grep "mapped (" | sed 's/.*(//g' | sed 's/ .*//g')
		duplicaterate=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiStats/$metrics | grep -A1 "LIBRARY"  | awk '{print $8}' | tail -1)
		coverage=$(echo $totalmapped | awk -v GenomeSize="$GenomeSize" '{print 100*$1/GenomeSize}')

		echo -e $Sample"\t"$totalreads"\t"$totalmapped"\t"$percmapped"\t"$duplicaterate"\t"$coverage >> ~/BM.Mapping.Summary.tsv
	done

	for j in $( \ls /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia | grep "bam" | grep "dedup");
	do
		Sample=$(echo $j | awk -F '.' '{print $1}')
		flagstat=$(echo $j | sed 's/bam/flagstat/g')
		metrics=$(echo $j | sed 's/dedup.bam/metrics/g')
		totalreads=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/$flagstat | head -1 | sed 's/ +.*//g')
		totalmapped=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/$flagstat | grep "mapped (" | sed 's/ +.*//g')
		percmapped=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/$flagstat | grep "mapped (" | sed 's/.*(//g' | sed 's/ .*//g')
		duplicaterate=$(cat /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/$metrics | grep -A1 "LIBRARY"  | awk '{print $9}' | tail -1)
		coverage=$(echo $totalmapped | awk -v GenomeSize="$GenomeSize" '{print 100*$1/GenomeSize}')

		echo -e $Sample"\t"$totalreads"\t"$totalmapped"\t"$percmapped"\t"$duplicaterate"\t"$coverage >> ~/BM.Mapping.Summary.tsv
	done

fi



if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/MergedVCF/")
	cd $outdir

	#PHASE THE DATA
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/beagle/beagle.jar gt=/local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.ref.vcf nthreads=32 ibd=true out=/local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle niterations=20

	##Fix sample IDs
	zcat /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.vcf.gz | sed 's/BM3/TRS-male_3-M1/g' | sed 's/BM4/TRS-male_4-M1/g' | sed 's/Bm_Malaysia/Bm-Malaysia/g' | sed 's/F3_male/F3-male/g' | sed 's/_male_M1/-male-M1/g'  | sed 's/_m1/-m1/g' | sed 's/_M1/-M1/g' | sed 's/TRS_male/TRS-male/g' | sed 's/W_male/W-male/g' | grep -E "Chr|#" > /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf
	/usr/local/packages/plink/plink2 --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf --make-king --homozyg --allow-extra-chr --out All.merged.phased.plink
	/usr/local/packages/plink/plink2 --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf --export 'ped' --allow-extra-chr --out All.merged.phased.plink
	/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf --recode --allow-extra-chr --out All.merged.phased.plink.1.9


	####Do Kinship on RoH regions only
	/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.RoH.vcf --allow-extra-chr --het --out All.merged.beagle.sample.RoH
	/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf --het --allow-extra-chr --out All.merged.phased.plink.1.9
	cat /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.ref.vcf | sed 's/BM3/TRS-male_3-M1/g' | sed 's/BM4/TRS-male_4-M1/g' | sed 's/Bm_Malaysia/Bm-Malaysia/g' | sed 's/F3_male/F3-male/g' | sed 's/_male_M1/-male-M1/g'  | sed 's/_m1/-m1/g' | sed 's/_M1/-M1/g' | sed 's/TRS_male/TRS-male/g' | sed 's/W_male/W-male/g' | grep -E "Chr|#" > /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.ref.header.vcf

	/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.ref.header.vcf --het --allow-extra-chr --out All.merged.chr

	/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.RoH.noX.vcf --allow-extra-chr --het --out All.merged.beagle.sample.RoH.noX
	cat All.merged.beagle.sample.RoH.het | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ~/All.merged.beagle.sample.RoH.het
fi


###Download the Brugia SRAs
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/BrugiaPopGenome_Fastqs/")
	rm -rf $outdir
	mkdir $outdir

	for i in $(cat ~/Brugia.malayi.populations.SRA.txt | awk '{print $2}' | sort | uniq)
	do
		mkdir /local/scratch/jmattick/BrugiaPopGenome_Fastqs/$i/
		cd /local/scratch/jmattick/BrugiaPopGenome_Fastqs/$i/
		for j in $(cat ~/Brugia.malayi.populations.SRA.txt | grep $i | awk '{print $1}' | sort | uniq)
		do

			qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/BM_PopGenome_Download.sh -i $i -j $j

		done
	done
fi

##Combine FASTQS
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BrugiaMalayiPopGenome_CombinedFastqs/")
	rm -rf $outdir
	mkdir $outdir
	#outdir=$(echo "/local/scratch/jmattick/BrugiaPopGenome_Fastqs/")

	for i in $(cat ~/Brugia.malayi.populations.SRA.txt | awk '{print $2}' | sort | uniq )
	do
		filename=$(echo $i | sed 's/\r//g')
		#mv $outdir/$i $outdir/$filename
		cat /local/scratch/jmattick/BrugiaPopGenome_Fastqs/$i/*_1.fastq > $outdir/$filename"_1.fastq"
		cat /local/scratch/jmattick/BrugiaPopGenome_Fastqs/$i/*_2.fastq > $outdir/$filename"_2.fastq"

		#rm -rf /local/scratch/jmattick/BrugiaPopGenome_Fastqs/$i/

	done
fi
##gzip FASTQS
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BrugiaMalayiPopGenome_CombinedFastqs/")

	for i in $(cat ~/Brugia.malayi.populations.SRA.txt | awk '{print $2}' | sort | uniq )
	do
		filename=$(echo $i | sed 's/\r//g')
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd gzip $outdir/$filename"_1.fastq"
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd gzip $outdir/$filename"_2.fastq"


	done

fi


###Download the Wuchereria SRAs
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/WBancroftiPopGenome_CombinedFastqs/")
	rm -rf $outdir
	mkdir $outdir

	for i in $(cat ~/Wuchereria.bancrofti.populations.SRA.txt | awk '{print $2}' | sort | uniq)
	do
		filename=$(echo $i | sed 's/\r//g')

		mkdir /local/scratch/jmattick/WuchPopGenome_Fastqs/$filename/
		cd /local/scratch/jmattick/WuchPopGenome_Fastqs/$filename/
		for j in $(cat ~/Wuchereria.bancrofti.populations.SRA.txt | grep $i | awk '{print $1}' | sort | uniq)
		do

			qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/WB_PopGenome_Download.sh -i $filename -j $j

		done
		filename=$(echo $i | sed 's/\r//g')
		#mv $outdir/$i $outdir/$filename

	done
fi


##Combine FASTQS for Wuch
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/WuchPopGenome_CombinedFastqs/")
	rm -rf $outdir
	mkdir $outdir

	for i in $(\ls /local/scratch/jmattick/WuchPopGenome_Fastqs/ | grep "")
	do
		filename=$(echo $i)
		#mv $outdir/$i $outdir/$filename
		#cat /local/scratch/jmattick/WuchPopGenome_Fastqs/$filename/*_1.fastq > $outdir/$filename"_1.fastq"
		#cat /local/scratch/jmattick/WuchPopGenome_Fastqs/$filename/*_2.fastq > $outdir/$filename"_2.fastq"

		#rm -rf /local/scratch/jmattick/WuchPopGenome_Fastqs/$filename/

		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/WuchCombine.sh -i $i

	done
	#qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/BMalayiClinical.sh 
fi

##Wuch GZIP
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/WuchPopGenome_CombinedFastqs/")

	for i in $(\ls /local/scratch/jmattick/WuchPopGenome_CombinedFastqs/ | grep "_1.fastq" )
	do
		fastq2=$(echo $i | sed 's/_1/_2/g')
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd gzip $outdir/$i
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd gzip $outdir/$fastq2


	done

fi

#Wuch Contig Mapping
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/BmtoWuchNucmer/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/mummer/nucmer -p NucmerBmtoWuch /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS14.genomic.fa
	/usr/local/packages/mummer/show-coords -qlTHb $outdir/NucmerBmtoWuch.delta > $outdir/NucmerBmtoWuch.coords


fi

##BMalayi Remap
if [[ 1 -eq 2 ]]
then

	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	#fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa)

	dictpath=$(echo $fastapath | sed 's/fa/dict/g')
	#samtools faidx $fastapath

	#/usr/local/packages/bwa/bin/bwa index $fastapath
	#rm $dictpath
	#/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath


	#outdir=$(echo "/local/scratch/jmattick/Malayi_PopGenomeRemap/")
	outdir=$(echo "/local/scratch/jmattick/Malayi_PopGenome_unmasked/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for fastq1 in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BrugiaMalayiPopGenome_CombinedFastqs/*_1.fastq.gz | grep "fastq");
	do
		j=$(echo $fastq1 | awk -F "/" '{print $NF}' | sed 's/_1.fastq//g')
		fastq2=$(echo $fastq1 | sed 's/_1.fastq.gz/_2.fastq.gz/g')

		qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

		list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
		
		list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

		list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools  depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf
		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --sample-ploidy 1 --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.haploid.vcf

	done
fi


if [[ 1 -eq 2 ]]
then

	outdir=$(echo "/local/scratch/jmattick/Malayi_Output/")

	rm -rf $outdir
	mkdir $outdir
	cp /local/scratch/jmattick/Malayi_PopGenomeRemap/*.sorted.dedup* $outdir/
	#cp /local/scratch/jmattick/Malayi_PopGenomeRemap/TRS_male_4.sorted.dedup* $outdir/





fi
##Wuch Fix
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/WuchPopGenome_CombinedPairedFastqs/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for fastq1 in $(\ls /local/scratch/jmattick/WuchPopGenome_CombinedFastqs/*_1.fastq.gz | grep "fastq");
	do
		fastq2=$(echo $fastq1 | sed 's/_1.fastq.gz/_2.fastq.gz/g')
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=30G -b y -cwd /home/jmattick/SNP_Scripts/fastqCombinePairedEnd.py $fastq1 $fastq2
	done
fi



##Wuch Map
if [[ 1 -eq 2 ]]
then

	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa)

	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	/usr/local/packages/bwa/bin/bwa index $fastapath
	samtools faidx $fastapath
	rm $dictpath
	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath


	outdir=$(echo "/local/scratch/jmattick/Wuch_PopGenomeRemap/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for fastq1 in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/WBancroftiPopGenome_CombinedFastqs/*_1.fastq.gz | grep "fastq");
	do
		j=$(echo $fastq1 | awk -F "/" '{print $NF}' | sed 's/_1.fastq.gz//g')
		fastq2=$(echo $fastq1 | sed 's/_1.fastq.gz/_2.fastq.gz/g')

		qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

		list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
		
		list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

		list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools  depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf

	done
fi


if [[ 1 -eq 2 ]]
then

	outdir=$(echo "/local/scratch/jmattick/Wuch_Output/")

	rm -rf $outdir
	mkdir $outdir
	cp /local/scratch/jmattick/Wuch_PopGenomeRemap/*.sorted.dedup* $outdir/

fi




##Final Copy
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/")
	rm -rf $outdir
	mkdir $outdir


	species="BMalayi"
	mkdir $outdir"/"$species

	cp /local/scratch/jmattick/FilteredVariants/$species/*.tsv $outdir"/"$species

	for i in $(\ls /local/scratch/jmattick/BMalayiMalaysia_Masked/ | grep "tsv")
	do
		out=$(echo $i | sed 's/_Masked//g')
		cp /local/scratch/jmattick/BMalayiMalaysia_Masked/$i $outdir"/"$species"/"$out
	done

	species="BPahangi"
	mkdir $outdir"/"$species

	cp /local/scratch/jmattick/FilteredVariants/$species/*.tsv $outdir"/"$species

	species="WBancrofti"
	mkdir $outdir"/"$species

	cp /local/scratch/jmattick/FilteredVariants/$species/*.tsv $outdir"/"$species

fi


#Final Residues

if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/")

	species="BMalayi"
	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	/home/jdhotopp/bin/residues.pl $fastapath > $outdir"/"$species".res.txt"

	species="BPahangi"
	fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)
	/home/jdhotopp/bin/residues.pl $fastapath > $outdir"/"$species".res.txt"


	species="WBancrofti"
	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS14.genomic_masked.fa)
	/home/jdhotopp/bin/residues.pl $fastapath > $outdir"/"$species".res.txt"

fi


if [[ 1 -eq 2 ]]
then
	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	outdir=$(echo "/local/scratch/jmattick/MalayiSNPs_Joint/")
	rm -rf $outdir
	mkdir $outdir
	filelist=$(ls /local/scratch/jmattick/*Malayi*/*.sorted.dedup.bam | tr '\n' ';' | sed 's/;/ --input /g' | sed 's/ --input $//g' | awk '{print "--input "$0}')
	j="Combined"
	cd $outdir
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath $filelist --output $outdir/$j.sorted.dedup.vcf

fi

###Download the OVolvulus SRAs
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OVolvulusPopGenome_CombinedFastqs/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for i in $(cat ~/OVolvulus.SRA.tsv | awk '{print $2}' | sort | uniq)
	do
		filename=$(cat ~/OVolvulus.SRA.tsv | grep $i | awk '{print $1}')

		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/OV_PopGenome_Download.sh -i $i -j $filename -outdir $outdir 

	done
fi

####Compress O Volvulus Fastqs
if [[ 1 -eq 2 ]]
then

	for i in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OVolvulusPopGenome_CombinedFastqs/*.fastq | grep ".fastq" )
	do
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd gzip $i


	done

fi

###Rename O Volvulus Samples
if [[ 1 -eq 2 ]]
then

	for i in $(cat ~/OVolvulus.SRA.tsv | awk '{print $2}' | sort | uniq )
	do
		filename=$(cat ~/OVolvulus.SRA.tsv | grep $i | awk '{print $1}')

		mv "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OVolvulusPopGenome_CombinedFastqs/"$i"_1.fastq.gz" "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OVolvulusPopGenome_CombinedFastqs/"$filename"_1.fastq.gz"
		mv "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OVolvulusPopGenome_CombinedFastqs/"$i"_2.fastq.gz" "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OVolvulusPopGenome_CombinedFastqs/"$filename"_2.fastq.gz"


	done

fi


###Create Records of BAM, depth, flagstat, metrics files
if [[ 1 -eq 2 ]]
then
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/B_malayi/
	cp /local/scratch/jmattick/*unmasked/*.sorted.dedup.b* /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/B_malayi/

fi

#######O Volvulus Mapping
if [[ 1 -eq 2 ]]
then
	#fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa)
	fastapath=$(\ls /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa)
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	/usr/local/packages/bwa/bin/bwa index $fastapath
	rm $dictpath
	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

	outdir=$(echo "/local/scratch/jmattick/OVolvulus_masked")
	#outdir=$(echo "/local/scratch/jmattick/BMalayiMalaysia")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for fastq1 in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OVolvulusPopGenome_CombinedFastqs/*_1.fastq.gz | grep "_1.fastq.gz");
	do
		j=$(echo $fastq1 | awk -F '/' '{print $NF}' | sed 's/_1.fastq.gz//g')
		fastq2=$(echo $fastq1 | sed 's/_1/_2/g')

		qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

		list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
		
		list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

		list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

		qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools  depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

		#qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk AddOrReplaceReadGroups --INPUT $outdir/$j.sorted.dedup.bam --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM $j --OUTPUT $outdir/$j.sorted.dedup.readgroup.bam --CREATE_INDEX true

		##Works. Why does BWA fuck it up?
		##samtools view -H $outdir/$j.sorted.dedup.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' |  samtools reheader - $outdir/$j.sorted.dedup.bam > $outdir/$j.sorted.dedup.reheader.bam

		#qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf
		#qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --sample-ploidy 1 --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.haploid.vcf

	done
fi


#######Joint Variant Calling directory
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/
	mkdir /local/scratch/jmattick/Joint_Genotyping/
fi


####Wuch Joint Genotyping Parallel freebayes
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/
	mkdir /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/
	cd /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/
	samtools faidx /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa
	\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/W_bancrofti/*.bam > /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/WB_filelist.txt
	/usr/local/packages/python-2.7.17/bin/python2.7 /local/scratch/jmattick/miniconda3/bin/fasta_generate_regions.py /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa 100000 > /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.regions
	qsub -V -P jhotopp-gcid-proj4b-filariasis -o /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/W_bancrofti_jointcall.vcf -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd freebayes-parallel /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.regions 36 -f /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa -L /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/WB_filelist.txt
	qsub -V -P jhotopp-gcid-proj4b-filariasis -o /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/W_bancrofti_jointcall.haploid.vcf -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd freebayes-parallel /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.regions 36 -f /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa -p 1 -L /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/WB_filelist.txt

fi


####B Malayi GVCF Calling
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/B_malayi/
	mkdir /local/scratch/jmattick/Joint_Genotyping/B_malayi/
	cd /local/scratch/jmattick/Joint_Genotyping/B_malayi/
	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	samtools faidx $fastapath

	for bam in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/B_malayi/*.bam | grep "bam");
	do
		samplename=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".g.vcf"}')
		samplenamehap=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".haploid.g.vcf"}')

		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --reference $fastapath --input $bam --output /local/scratch/jmattick/Joint_Genotyping/B_malayi/$samplename
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --sample-ploidy 1 --reference $fastapath --input $bam --output /local/scratch/jmattick/Joint_Genotyping/B_malayi/$samplenamehap
	done
fi

####B Pahangi GVCF Calling
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/B_pahangi/
	mkdir /local/scratch/jmattick/Joint_Genotyping/B_pahangi/
	cd /local/scratch/jmattick/Joint_Genotyping/B_pahangi/
	fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)
	dictpath=$(echo $fastapath | sed 's/fasta/dict/g')
	samtools faidx $fastapath
	#/usr/local/packages/bwa/bin/bwa index $fastapath
	rm $dictpath
	/usr/local/packages/jdk-8u171/bin/java -Xmx12g -jar /usr/local/packages/picard/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

	for bam in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/B_pahangi/*.bam | grep "bam" | grep -v "\-FR3");
	do
		samplename=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".g.vcf"}')
		samplenamehap=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".haploid.g.vcf"}')

		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --reference $fastapath --input $bam --output /local/scratch/jmattick/Joint_Genotyping/B_pahangi/$samplename
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --sample-ploidy 1 --reference $fastapath --input $bam --output /local/scratch/jmattick/Joint_Genotyping/B_pahangi/$samplenamehap
	done
fi

####O Volvulus GVCF Calling
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/
	mkdir /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/
	cd /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/
	fastapath=$(\ls /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa)
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')
	samtools faidx $fastapath
	/usr/local/packages/bwa/bin/bwa index $fastapath
	rm $dictpath
	/usr/local/packages/jdk-8u171/bin/java -Xmx12g -jar /usr/local/packages/picard/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

	for bam in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/O_Volvulus/*.bam | grep "bam" );
	do
		samplename=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".g.vcf"}')
		samplenamehap=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".haploid.g.vcf"}')

		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --reference $fastapath --input $bam --output /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/$samplename
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --sample-ploidy 1 --reference $fastapath --input $bam --output /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/$samplenamehap
	done
fi


if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/MummerPlots/
	mkdir /local/scratch/jmattick/Mummerplots/
	cd /local/scratch/jmattick/Mummerplots/

	MalayiPath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa)
	/usr/local/packages/mummer/nucmer -p BMvsWB $MalayiPath $fastapath
	/usr/local/packages/mummer/show-coords -qlTH BMvsWB.delta > BMvsWB.coords

	fastapath=$(\ls /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa)
	/usr/local/packages/mummer/nucmer -p BMvsOV $MalayiPath $fastapath
	/usr/local/packages/mummer/show-coords -qlTH BMvsOV.delta > BMvsOV.coords


fi

####OV GVCF DB Construction
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/DB/
	cd /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/
	IntervalList=$(cat /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

	GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/*.g.vcf | grep -v "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/DB/ $IntervalList

	rm -rf /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/DB_haploid/
	cd /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/
	IntervalList=$(cat /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

	GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/*.g.vcf | grep "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/DB_haploid/ $IntervalList

fi

####BM GVCF DB Construction
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/B_malayi/DB/
	cd /local/scratch/jmattick/Joint_Genotyping/B_malayi/
	IntervalList=$(cat /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

	GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/B_malayi/*.g.vcf | grep -v "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/B_malayi/DB/ $IntervalList

	rm -rf /local/scratch/jmattick/Joint_Genotyping/B_malayi/DB_haploid/
	cd /local/scratch/jmattick/Joint_Genotyping/B_malayi/
	IntervalList=$(cat /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

	GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/B_malayi/*.g.vcf | grep "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/B_malayi/DB_haploid/ $IntervalList

fi

####BP GVCF DB Construction
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/Joint_Genotyping/B_pahangi/DB/
	cd /local/scratch/jmattick/Joint_Genotyping/B_pahangi/
	IntervalList=$(cat /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

	GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/B_pahangi/*.g.vcf | grep -v "BEI_FR3_Bpahangi_multi_AF_a" | grep -v "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
	list=$(qstat | grep "gatk" | grep -v "399411" | grep -v "hqw" | awk '{print $1}' | paste -s -d,)

	qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/B_pahangi/DB/ $IntervalList

	rm -rf /local/scratch/jmattick/Joint_Genotyping/B_pahangi/DB_haploid/
	cd /local/scratch/jmattick/Joint_Genotyping/B_pahangi/
	IntervalList=$(cat /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

	GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/B_pahangi/*.g.vcf | grep -v "BEI_FR3_Bpahangi_multi_AF_a" | grep "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
	list=$(qstat | grep "gatk" | grep -v "399411" | grep -v "hqw" | awk '{print $1}' | paste -s -d,)

	qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/B_pahangi/DB_haploid/ $IntervalList

fi


###Variant Calling from DBs
if [[ 1 -eq 2 ]]
then
	##OV
	fastapath=$(\ls /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa)
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/O_Volvulus/DB/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/OV.final.jointcall.vcf 
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/O_Volvulus/DB_haploid/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/OV.final.jointcall.haploid.vcf 
	##BM
	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/B_malayi/DB/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/B_malayi/BM.final.jointcall.vcf 
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/B_malayi/DB_haploid/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/B_malayi/BM.final.jointcall.haploid.vcf 
	##BP
	fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/B_pahangi/DB/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/B_pahangi/BP.final.jointcall.vcf 
	qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/B_pahangi/DB_haploid/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/B_pahangi/BP.final.jointcall.haploid.vcf 

fi

if [[ 1 -eq 2 ]]
then
	~sdaugherty/bin/reference_genome_match_t iler.pl --mummer_coords_file amos_newbler.coords --min_match_length 100 --mummer_delta_file amos_newbler.delta --method nucmer > result2.delta
	mummerplot -Qfile ordered_query_list.txt result2.delta
	cat ordered_query_list.txt | perl -ne 'chomp; ($asmbl,$len,$dir) = split/\t/; print "EWANA.pseudomolecule.3\t$asmbl\t$dir\n";' > test2
	~sdaugherty/svn/jorvis_utilities/general/create_fasta_pseudomolecules.pl --input_fasta_file=454LargeContigs_over3kbp.fna --map_file=test2 --output_file=EWANA.pseudomolecule.3.fasta --log=log.txt --unmapped_output=unmapped.fsa
	#Didn't work...
fi

if [[ 1 -eq 2 ]]
then
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/W_bancrofti_jointcall.vcf --window-pi 10000 --out /local/scratch/jmattick/W_bancrofti_jointcall.pi
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/OV.final.jointcall.vcf --window-pi 10000 --out /local/scratch/jmattick/OV.final.jointcall.pi
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/B_malayi/BM.final.jointcall.vcf --window-pi 10000 --out /local/scratch/jmattick/BM.final.jointcall.pi
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/B_pahangi/BP.final.jointcall.vcf --window-pi 10000 --out /local/scratch/jmattick/BP.final.jointcall.pi


fi

if [[ 1 -eq 2 ]]
then
	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa)
	inputvcf=$(\ls /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/W_bancrofti_jointcall.vcf)
	output=$(echo $inputvcf | sed 's/vcf/filtered.vcf/g')
	outputfiltered=$(echo $inputvcf | sed 's/vcf/filteredonly.vcf/g')
	#/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"
	/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-name "QD" --filter-expression "QD < 5.0" --filter-name "QUAL" --filter-expression "QUAL < 30.0" --filter-name "DP" --filter-expression "DP < 14.0" --filter-name "MQ" --filter-expression "MQ < 30.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "FS" --filter-expression "FS > 60.0" 

	/usr/local/packages/gatk-4.1.7.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outputfiltered -select 'vc.isNotFiltered()'

	fastapath=$(\ls /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa)
	inputvcf=$(\ls /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/OV.final.jointcall.vcf)
	output=$(echo $inputvcf | sed 's/vcf/filtered.vcf/g')
	outputfiltered=$(echo $inputvcf | sed 's/vcf/filteredonly.vcf/g')
	#/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"
	/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-name "QD" --filter-expression "QD < 5.0" --filter-name "QUAL" --filter-expression "QUAL < 30.0" --filter-name "DP" --filter-expression "DP < 14.0" --filter-name "MQ" --filter-expression "MQ < 30.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "FS" --filter-expression "FS > 60.0" 

	/usr/local/packages/gatk-4.1.7.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outputfiltered -select 'vc.isNotFiltered()'


	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	inputvcf=$(\ls /local/scratch/jmattick/Joint_Genotyping/B_malayi/BM.final.jointcall.vcf)
	output=$(echo $inputvcf | sed 's/vcf/filtered.vcf/g')
	outputfiltered=$(echo $inputvcf | sed 's/vcf/filteredonly.vcf/g')
	#/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"
	/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-name "QD" --filter-expression "QD < 5.0" --filter-name "QUAL" --filter-expression "QUAL < 30.0" --filter-name "DP" --filter-expression "DP < 14.0" --filter-name "MQ" --filter-expression "MQ < 30.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "FS" --filter-expression "FS > 60.0" 

	/usr/local/packages/gatk-4.1.7.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outputfiltered -select 'vc.isNotFiltered()'

	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	inputvcf=$(\ls /local/scratch/jmattick/Joint_Genotyping/B_malayi/BM.final.jointcall.haploid.vcf)
	output=$(echo $inputvcf | sed 's/vcf/filtered.vcf/g')
	outputfiltered=$(echo $inputvcf | sed 's/vcf/filteredonly.vcf/g')
	#/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"
	/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-name "QD" --filter-expression "QD < 5.0" --filter-name "QUAL" --filter-expression "QUAL < 30.0" --filter-name "DP" --filter-expression "DP < 7.0" --filter-name "MQ" --filter-expression "MQ < 30.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "FS" --filter-expression "FS > 60.0" 

	/usr/local/packages/gatk-4.1.7.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outputfiltered -select 'vc.isNotFiltered()'



	fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)
	inputvcf=$(\ls /local/scratch/jmattick/Joint_Genotyping/B_pahangi/BP.final.jointcall.vcf )
	output=$(echo $inputvcf | sed 's/vcf/filtered.vcf/g')
	outputfiltered=$(echo $inputvcf | sed 's/vcf/filteredonly.vcf/g')
	#/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"
	/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-name "QD" --filter-expression "QD < 5.0" --filter-name "QUAL" --filter-expression "QUAL < 30.0" --filter-name "DP" --filter-expression "DP < 14.0" --filter-name "MQ" --filter-expression "MQ < 30.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "FS" --filter-expression "FS > 60.0" 
	/usr/local/packages/gatk-4.1.7.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outputfiltered -select 'vc.isNotFiltered()'

fi


###Split B. malayi samples back out
if [[ 1 -eq 2 ]]
then

fi

if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/PopGenDepthFiles/")
	rm -rf $outdir
	mkdir $outdir

	for bam in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/*/*.bam)
	do
		sample=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		species=$(echo $bam | awk -F '/' '{print $7}' | awk -F '.' '{print $1}')
		qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$species.$sample.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools depth -aa -d 10000000000 $bam

	done

	for bam in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/All_Nematode_BAMs/*/*.bam | grep -v "W_ban")
	do
		sample=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		species=$(echo $bam | awk -F '/' '{print $7}' | awk -F '.' '{print $1}')
		qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$species.$sample.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools depth -aa -d 10000000000 $bam

	done

fi

if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenDepthFile_table/")
	rm -rf $outdir
	mkdir $outdir

	cd $outdir

	for depth in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenDepthFiles/*.depth | grep "depth")
	do
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/DepthPlotCalculator_mod.sh -depth $depth -outdir $outdir 
	done
fi


if [[ 1 -eq 2 ]]
then
	#for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep -v "B_" | grep -v "W_" | grep -v "Diro" | grep -v "O_" | grep -v "mixed")
	for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep -v "B_" | grep -v "W_" | grep "Wban" | grep -v "O_" | grep -v "mixed")
	do
		fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')

		inputvcf=$(\ls /local/scratch/jmattick/Joint_Genotyping/$species/$species.final.jointcall.vcf)
		output=$(echo $inputvcf | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $inputvcf | sed 's/vcf/filteredonly.vcf/g')
		#/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"
		/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-name "QD" --filter-expression "QD < 5.0" --filter-name "QUAL" --filter-expression "QUAL < 30.0" --filter-name "DP" --filter-expression "DP < 14.0" --filter-name "MQ" --filter-expression "MQ < 30.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "FS" --filter-expression "FS > 60.0" 
   
		#/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"
		/usr/local/packages/gatk-4.1.7.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $output --filter-name "QD" --filter-expression "QD < 5.0" --filter-name "QUAL" --filter-expression "QUAL < 30.0" --filter-name "DP" --filter-expression "DP < 14.0" --filter-name "MQ" --filter-expression "MQ < 30.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "FS" --filter-expression "FS > 60.0" 
   
		/usr/local/packages/gatk-4.1.7.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outputfiltered -select 'vc.isNotFiltered()'
	done

fi
if [[ 1 -eq 2 ]]
then
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/W_bancrofti_jointcall.filteredonly.vcf --window-pi 10000 --out /local/scratch/jmattick/W_bancrofti_jointcall.filteredonly.pi
	#Excluding Samples
	/usr/local/packages/vcftools/bin/vcftools --indv Haiti1007-4 --indv Haiti1814-4 --indv Haiti2070-1 --indv Kenya0258-1 --indv Mali0159-12 --vcf /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/W_bancrofti_jointcall.filteredonly.vcf --window-pi 10000 --out /local/scratch/jmattick/W_bancrofti_jointcall.filteredonly.subset.pi

	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/O_Volvulus/OV.final.jointcall.filteredonly.vcf --window-pi 10000 --out /local/scratch/jmattick/OV.final.jointcall.filteredonly.pi
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/B_malayi/BM.final.jointcall.filteredonly.vcf --window-pi 10000 --out /local/scratch/jmattick/BM.final.jointcall.filteredonly.pi
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/B_pahangi/BP.final.jointcall.filteredonly.vcf --window-pi 10000 --out /local/scratch/jmattick/BP.final.jointcall.filteredonly.pi
	cp /local/scratch/jmattick/*.filteredonly.* /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/
fi

if [[ 1 -eq 2 ]]
then
	#for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep -v "B_" | grep -v "W_" | grep -v "Wb" | grep -v "O_" | grep -v "mixed")
	for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep -v "B_" | grep -v "W_" | grep "Wban" | grep -v "O_" | grep -v "mixed")
	do
		fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')

		inputvcf=$(\ls /local/scratch/jmattick/Joint_Genotyping/$species/$species.final.jointcall.vcf)
		output=$(echo $inputvcf | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $inputvcf | sed 's/vcf/filteredonly.vcf/g')

		/usr/local/packages/vcftools/bin/vcftools --vcf $outputfiltered --window-pi 10000 --out /local/scratch/jmattick/$species.final.jointcall.filteredonly.pi

	done	
fi

##Bi-allelic WB
if [[ 1 -eq 2 ]]
then
	/usr/local/packages/vcftools/bin/vcftools --vcf /local/scratch/jmattick/Joint_Genotyping/W_bancrofti_freebayes/W_bancrofti_jointcall.filteredonly.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out /local/scratch/jmattick/W_bancrofti_jointcall.filteredonly.biallelic
	cat /local/scratch/jmattick/W_bancrofti_jointcall.filteredonly.biallelic.recode.vcf | grep -v "#" | awk '{print $1"\t"$2"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25}' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' | sed 's/\.:\.:\.:\.:\.:\.:\./2/g' > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/W_bancrofti/W_bancrofti_jointcall.filteredonly.biallelic.recode.vcfprocessed
fi



####BIG species download
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/")
	#rm -rf $outdir
	#mkdir $outdir
	#for sralist in $(\ls ~/*.sralist.txt)
	for sralist in $(\ls ~/Wbancrofti.highdepth.sra.txt | grep "Wban")
	do
		species=$(echo $sralist | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		rm -rf $outdir/$species
		mkdir $outdir/$species
		for sample in $(cat $sralist | awk '{print $2}' | sed 's/\r//g' | sort | uniq)
		do

			mkdir $outdir/$species/$sample/
			cd $outdir/$species/$sample/
			for sra in $(cat $sralist | grep $sample | awk '{print $1}')
			do
				qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/General_PopGenome_Download.sh -sra $sra -species $species -sample $sample -outdir $outdir 
			done
		done

	done

fi

##Big Compression
if [[ 1 -eq 2 ]]
then
	cd /local/scratch/jmattick/
	for i in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/*/*/*.fastq | grep ".fastq" )
	do
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd gzip $i


	done

fi


##Big Mapping
if [[ 1 -eq 2 ]]
then
	#rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/
	#mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/
	rm -rf /local/scratch/jmattick/temp/
	for species in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/ | grep "Wban" )
	do
		fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')
		dictpath=$(echo $fastapath | sed 's/fa/dict/g')

		#/usr/local/packages/bwa/bin/bwa index $fastapath
		#rm $dictpath
		#/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

		outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/"$species)
		#outdir=$(echo "/local/scratch/jmattick/BMalayiMalaysia")
		rm -rf $outdir
		mkdir $outdir
		cd $outdir

		for sample in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species | grep "");
		do
			fastqnumber=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_1.fastq.gz | wc -l)
			if [[ $fastqnumber -gt 1 ]]
			then
				cat /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_1.fastq.gz > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_1.fastq.gz
				cat /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_2.fastq.gz > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_2.fastq.gz
				fastq1=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_1.fastq.gz)
				fastq2=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_2.fastq.gz)
			else
				fastq1=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_1.fastq.gz | grep "fastq")
				fastq2=$(echo $fastq1 | sed 's/_1.fastq.gz/_2.fastq.gz/g')
			fi
			j=$sample

			qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

			list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
			
			list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt
			#rm $outdir/$j.sorted.metrics
			#rm $outdir/$j.sorted.dedup.bam
			#rm $outdir/$j.sorted.dedup.bai
			#qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=30G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=1000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 > $outdir/$j.BPDedup.txt

		done
	done


fi


if [[ 1 -eq 2 ]]
then
	for i in $(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/*/*.ba* | grep -v "dedup")
	do 
		rm $i
	done
	#mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/Celegans
	#cp /local/scratch/jmattick/CElegans_temp/* /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/Celegans/
	#rm -rf /local/scratch/jmattick/CElegans_temp/
fi



##Variant Calling for all except Diro

if [[ 1 -eq 2 ]]
then
	#for species in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/ | grep -v "Diro" )
	for species in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/ )
	do
		fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')
		dictpath=$(echo $fastapath | sed 's/fa/dict/g')
		samtools faidx $fastapath
		/usr/local/packages/bwa/bin/bwa index $fastapath
		rm $dictpath
		/usr/local/packages/jdk-8u171/bin/java -Xmx12g -jar /usr/local/packages/picard/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

		outdir=$(echo "/local/scratch/jmattick/Joint_Genotyping/"$species)
		rm -rf $outdir
		mkdir $outdir
		cd $outdir
		for bam in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/*.bam | grep "bam" );
		do
			samplename=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".g.vcf"}')
			samplenamehap=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".haploid.g.vcf"}')

			qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --reference $fastapath --input $bam --output $outdir/$samplename
			qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --sample-ploidy 1 --reference $fastapath --input $bam --output $outdir/$samplenamehap
		done
	done

fi

##Variant Calling for Diro -- Di is blood (mixed), JS is adults from Sydney
if [[ 1 -eq 2 ]]
then
	species=$(echo "Dirofilariaimmitis")
	fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')
	samtools faidx $fastapath
	#/usr/local/packages/bwa/bin/bwa index $fastapath
	#rm $dictpath
	#/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

	outdir=$(echo "/local/scratch/jmattick/Joint_Genotyping/"$species)
	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	for bam in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/*.bam | grep "bam" | grep "JS");
	do
		samplename=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".g.vcf"}')
		samplenamehap=$(echo $bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1".haploid.g.vcf"}')

		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --reference $fastapath --input $bam --output $outdir/$samplename
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --emit-ref-confidence GVCF --sample-ploidy 1 --reference $fastapath --input $bam --output $outdir/$samplenamehap
	done

	outdir=$(echo "/local/scratch/jmattick/Joint_Genotyping/"$species"_mixedsamples")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	regionpath=$(echo $fastapath | awk -F '/' '{print $NF}' | sed 's/fa/regions/g')

	\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/Di*.bam > $outdir/$species.filelist.txt
	/usr/local/packages/python-2.7.14/bin/python2.7 /local/scratch/jmattick/miniconda3/bin/fasta_generate_regions.py $fastapath 100000 > $outdir/$regionpath
	qsub -V -P jhotopp-gcid-proj4b-filariasis -o $outdir/Mixed_jointcall.vcf -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd freebayes-parallel $outdir/$regionpath 36 -f $fastapath -L $outdir/$species.filelist.txt
	qsub -V -P jhotopp-gcid-proj4b-filariasis -o $outdir/Mixed_jointcall.haploid.vcf -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd freebayes-parallel $outdir/$regionpath 36 -f $fastapath -p 1 -L $outdir/$species.filelist.txt


fi


####Nigon Element Determination: WB
if [[ 1 -eq 2 ]]
then
	onchofastapath=$(\ls /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa)
	brugiafastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	celegansfastapath=$(\ls /local/scratch/jmattick/caenorhabditis_elegans.PRJNA13758.WBPS15.genomic_masked.fa)
	rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/
	

	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/W_bancrofti/
	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa)
	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/W_bancrofti/
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p W_bancrofti_vs_BM $fastapath $brugiafastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p W_bancrofti_vs_OV $fastapath $onchofastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p W_bancrofti_vs_CE $fastapath $celegansfastapath




	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/B_pahangi/
	fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)
	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/B_pahangi/
  	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p B_pahangi_vs_BM $fastapath $brugiafastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p B_pahangi_vs_OV $fastapath $onchofastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p B_pahangi_vs_CE $fastapath $celegansfastapath

fi

if [[ 1 -eq 2 ]]
then
	onchofastapath=$(\ls /local/scratch/jmattick/onchocerca_volvulus.PRJEB513.WBPS15.genomic_masked.fa)
	brugiafastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)
	celegansfastapath=$(\ls /local/scratch/jmattick/caenorhabditis_elegans.PRJNA13758.WBPS15.genomic_masked.fa)

	rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/L_loa/
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/L_loa/
	fastapath=$(\ls /local/scratch/jmattick/loa_loa.PRJNA246086.WBPS15.genomic_masked.fa)
	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/L_loa/
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p L_loa_vs_BM $fastapath $brugiafastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p L_loa_vs_OV $fastapath $onchofastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p L_loa_vs_CE $fastapath $celegansfastapath



	rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/D_immitis/
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/D_immitis/
	fastapath=$(\ls /local/scratch/jmattick/dirofilaria_immitis.PRJEB1797.WBPS15.genomic_masked.fa)
	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/D_immitis/
  	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p D_immitis_vs_BM $fastapath $brugiafastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p D_immitis_vs_OV $fastapath $onchofastapath
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd nucmer -p D_immitis_vs_CE $fastapath $celegansfastapath

fi


if [[ 1 -eq 2 ]]
then
	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/W_bancrofti/
	/usr/local/packages/mummer/show-coords -qlTHb W_bancrofti_vs_BM.delta > W_bancrofti_vs_BM.coords
	/usr/local/packages/mummer/show-coords -qlTHb W_bancrofti_vs_OV.delta > W_bancrofti_vs_OV.coords
	/usr/local/packages/mummer/show-coords -qlTHb W_bancrofti_vs_CE.delta > W_bancrofti_vs_CE.coords
	cp /local/scratch/jmattick/W_bancrofti_jointcall.filteredonly.pi.windowed.pi ./
	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS15.genomic_masked.fa)
	/home/jdhotopp/bin/residues.pl $fastapath > W_bancrofti.res.txt

	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/B_pahangi/
	/usr/local/packages/mummer/show-coords -qlTHb B_pahangi_vs_BM.delta > B_pahangi_vs_BM.coords
	/usr/local/packages/mummer/show-coords -qlTHb B_pahangi_vs_OV.delta > B_pahangi_vs_OV.coords
	/usr/local/packages/mummer/show-coords -qlTHb B_pahangi_vs_CE.delta > B_pahangi_vs_CE.coords
	cp /local/scratch/jmattick/BP.final.jointcall.filteredonly.pi.windowed.pi ./
	fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)
	/home/jdhotopp/bin/residues.pl $fastapath > B_pahangi.res.txt


fi

if [[ 1 -eq 2 ]]
then
	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/LoaLoa/
	/usr/local/packages/mummer/show-coords -qlTHb L_loa_vs_BM.delta > LoaLoa_vs_BM.coords
	/usr/local/packages/mummer/show-coords -qlTHb L_loa_vs_OV.delta > Loaloa_vs_OV.coords
	/usr/local/packages/mummer/show-coords -qlTHb L_loa_vs_CE.delta > LoaLoa_vs_CE.coords
	cp /local/scratch/jmattick/LoaLoa.final.jointcall.filteredonly.pi.windowed.pi ./
	fastapath=$(\ls /local/scratch/jmattick/loa_loa.PRJNA246086.WBPS15.genomic_masked.fa)
	/home/jdhotopp/bin/residues.pl $fastapath > LoaLoa.res.txt

	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/Dirofilariaimmitis/
	/usr/local/packages/mummer/show-coords -qlTHb D_immitis_vs_BM.delta > Dirofilariaimmitis_vs_BM.coords
	/usr/local/packages/mummer/show-coords -qlTHb D_immitis_vs_OV.delta > Dirofilariaimmitis_vs_OV.coords
	/usr/local/packages/mummer/show-coords -qlTHb D_immitis_vs_CE.delta > Dirofilariaimmitis_vs_CE.coords
	cp /local/scratch/jmattick/Dirofilariaimmitis.final.jointcall.filteredonly.pi.windowed.pi ./
	fastapath=$(\ls /local/scratch/jmattick/dirofilaria_immitis.PRJEB1797.WBPS15.genomic_masked.fa)
	/home/jdhotopp/bin/residues.pl $fastapath > Dirofilariaimmitis.res.txt

	rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/Celegans
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/Celegans

	rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/Dmelanogaster
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/Dmelanogaster

	cp /local/scratch/jmattick/Celegans.final.jointcall.filteredonly.pi.windowed.pi /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/Celegans
	cp /local/scratch/jmattick/Dmelanogaster.final.jointcall.filteredonly.pi.windowed.pi /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/Dmelanogaster
fi


####General DB Construction
if [[ 1 -eq 2 ]]
then
	list=$(qstat | grep "gatk" | awk '{print $1}' | tr '\n' ',' | sed 's/,$//g')
	#for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep -v "B_" | grep -v "W_" | grep -v "mixed" | grep -v "O_")
	rm -rf /local/scratch/jmattick/gatktemp/
	mkdir /local/scratch/jmattick/gatktemp/
	for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep "B_" | grep -v "mixed")
	do
		rm -rf /local/scratch/jmattick/Joint_Genotyping/$species/DB/
		cd /local/scratch/jmattick/Joint_Genotyping/$species/
		fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')

		IntervalList=$(cat $fastapath | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

		GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/$species/*.g.vcf | grep -v "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
		qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/$species/DB/ $IntervalList --tmp-dir /local/scratch/jmattick/gatktemp/
		#qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/$species/DB/ $IntervalList --tmp-dir /local/scratch/jmattick/gatktemp/

		rm -rf /local/scratch/jmattick/Joint_Genotyping/$species/DB_haploid/
		cd /local/scratch/jmattick/Joint_Genotyping/$species/
		IntervalList=$(cat $fastapath | grep ">" | awk '{print $1}' | sed 's/>//g' | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -L /g' | awk '{print "-L "$0}')

		GVCFList=$(\ls /local/scratch/jmattick/Joint_Genotyping/$species/*.g.vcf | grep "haploid" | tr '\n' ' ' | sed 's/ $//g' | sed 's/ / -V /g' | awk '{print "-V "$0}')
		qsub -hold_jid $list -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/$species/DB_haploid/ $IntervalList --tmp-dir /local/scratch/jmattick/gatktemp/
		#qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenomicsDBImport $GVCFList --genomicsdb-workspace-path /local/scratch/jmattick/Joint_Genotyping/$species/DB_haploid/ $IntervalList --tmp-dir /local/scratch/jmattick/gatktemp/

	done
fi

if [[ 1 -eq 2 ]]
then
	##All species variant calling
	#for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep -v "B_" | grep -v "W_" | grep -v "Diro" | grep -v "O_" | grep -v "mixed")
	for species in $(\ls /local/scratch/jmattick/Joint_Genotyping/ | grep "Wban" | grep -v "mixed")
	do
		fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/$species/DB/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/$species/$species.final.jointcall.vcf 
		qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.1.7.0/gatk GenotypeGVCFs -R $fastapath -V gendb:///local/scratch/jmattick/Joint_Genotyping/$species/DB_haploid/ --allow-old-rms-mapping-quality-annotation-data -O /local/scratch/jmattick/Joint_Genotyping/$species/$species.final.jointcall.haploid.vcf 
	done
fi


###Mansonella mapping
##Big Mapping
if [[ 1 -eq 2 ]]
then
	#rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/
	#mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/
	rm -rf /local/scratch/jmattick/temp/
	species="Mperstans"
	fastapath=$(cat ~/PopGenome.fasta.specieslist.txt | grep $species | awk '{print $2}')
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	/usr/local/packages/bwa/bin/bwa index $fastapath
	rm $dictpath
	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/"$species)
	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for sample in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species | grep "");
	do
		fastqnumber=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_1.fastq.gz | wc -l)
		if [[ $fastqnumber -gt 1 ]]
		then
			cat /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_1.fastq.gz > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_1.fastq.gz
			cat /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_2.fastq.gz > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_2.fastq.gz
			fastq1=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_1.fastq.gz)
			fastq2=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome_Mapping/$species/$sample.combined_2.fastq.gz)
		else
			fastq1=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OtherSpeciesPopGenome/$species/$sample/*_1.fastq.gz | grep "fastq")
			fastq2=$(echo $fastq1 | sed 's/_1.fastq.gz/_2.fastq.gz/g')
		fi
		j=$sample

			#qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

			#list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			#qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
			
			#list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			#qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt
		rm $outdir/$j.sorted.metrics
		rm $outdir/$j.sorted.dedup.bam
		rm $outdir/$j.sorted.dedup.bai
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=30G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=1000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 > $outdir/$j.BPDedup.txt

	done

fi


