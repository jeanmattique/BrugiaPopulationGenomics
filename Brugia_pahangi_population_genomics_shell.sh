
##Working Directory Setup
if [[ 1 -eq 2 ]]
then

	rm -rf /local/scratch/jmattick/PahangiIllumina/
	mkdir /local/scratch/jmattick/PahangiIllumina/
fi


if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/PahangiIllumina_Masked/
	mkdir /local/scratch/jmattick/PahangiIllumina_Masked/
fi


if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/temp/
	mkdir /local/scratch/jmattick/temp/
fi


####FASTA File Masking
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/Pahangi_FinalRepeat/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	cp /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.fasta $outdir

	/usr/local/packages/repeatmodeler-1.0.11/BuildDatabase -name BrugiaPahangi.FINAL.V5.2.fasta -engine wublast $outdir/BrugiaPahangi.FINAL.V5.2.fasta

	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/repeatmodeler-1.0.11/RepeatModeler -database BrugiaPahangi.FINAL.V5.2.fasta > $outdir/RepeatModeler.txt

	list=$(cat $outdir/RepeatModeler.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/repeatmasker-4.0.7/RepeatMasker -lib $outdir/BrugiaPahangi.FINAL.V5.2.fasta-families.fa $outdir/BrugiaPahangi.FINAL.V5.2.fasta


fi
#fastapath=$(\ls /local/aberdeen2rw/julie/JM_dir/PahangiPilonFASTA/1.BrugiaPahangi.CURRENT.fasta)
#fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.fasta)
fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)

dictpath=$(echo $fastapath | sed 's/fasta/dict/g')
mmipath=$(echo $fastapath | sed 's/fasta/mmi/g')


if [[ 1 -eq 2 ]]
then
	rm -rf $dictpath
	/usr/local/packages/bwa/bin/bwa index $fastapath
	/usr/local/packages/samtools/bin/samtools faidx $fastapath
	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath
	minimap2 -d $mmipath $fastapath
fi

if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/MultiAF/")
	#outdir=$(echo "/local/scratch/jmattick/PahangiIllumina_Masked/MultiAF/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	##Location of multi adult female B pahangi Sample: /local/projects/EBMAL/
	for j in $( \ls /local/projects/EBMAL/ | grep "BEI_FR3_Bpahangi_multi_AF");
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

	done

fi


##Download and Set up FastQs
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ")
	rm -rf $outdir
	mkdir $outdir
	for i in $(cat /home/jmattick/BPahangiSRA.tsv | awk '{print $1}' | sort | uniq)
	do
		SRANumber=$(cat /home/jmattick/BPahangiSRA.tsv | grep $i | wc -l)
		if [[ $SRANumber -eq 2 ]]
		then
			SRA1=$(cat /home/jmattick/BPahangiSRA.tsv | grep $i | awk '{print $2}' | head -1)
			SRA2=$(cat /home/jmattick/BPahangiSRA.tsv | grep $i | awk '{print $2}' | tail -1)

			qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /home/jmattick/SNP_Scripts/Pahangi_FASTQ_Combine.sh -SRA1 $SRA1 -SRA2 $SRA2 -Sample $i
		else
			SRA1=$(cat /home/jmattick/BPahangiSRA.tsv | grep $i | awk '{print $2}')
			qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /home/jmattick/SNP_Scripts/Pahangi_FASTQ_Combine.sh -SRA1 $SRA1 -SRA2 "none" -Sample $i

		fi

	done
fi


if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/Clinical/")
	#outdir=$(echo "/local/scratch/jmattick/PahangiIllumina_Masked/Clinical/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	##Location of clinical B pahangi Samples: /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/


	for fastq1 in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/*_1.fastq.gz | grep "clinical");
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
	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/FR3_old/")
	#outdir=$(echo "/local/scratch/jmattick/PahangiIllumina_Masked/FR3_old/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	##Location of first FR3 B pahangi Samples: /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/

	for fastq1 in $( \ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/*_1.fastq.gz | grep "FR3");
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
	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/FR3_UWO/")
	#outdir=$(echo "/local/scratch/jmattick/PahangiIllumina_Masked/FR3_UWO/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	##Location of second FR3 B pahangi Samples: /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/

	for j in $( \ls /local/projects/EBMAL/ | grep "FR3_UWO");
	do
		fastqnumber=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | wc -l )
		if [[ $fastqnumber -eq 2 ]]
		then
			fastq1=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | head -1)
			fastq2=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | tail -1)

			qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /home/jmattick/SNP_Scripts/Pahangi_FASTQ_Combine_EBMAL.sh -fastq1 $fastq1 -fastq2 $fastq2 -Sample $j > $outdir/$j.Combine.txt

			list=$(cat $outdir/$j.Combine.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			output=$(echo $j"_R1.fastq.gz")

			fastq1=$(echo /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output)
			output=$(echo $j"_R2.fastq.gz")

			fastq2=$(echo /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output)


			#qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

			list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
			
			list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

			list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf


		else
			fastq1=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | grep "fastq")
			fastq2=$(echo $fastq1 | sed 's/R1.fastq.gz/R2.fastq.gz/g')


			qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

			list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
			
			list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

			list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf

		fi
	done

fi


if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/Bp1AM/")
	#outdir=$(echo "/local/scratch/jmattick/PahangiIllumina_Masked/Bp1AM/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	##Location of final FR3 B pahangi Samples: /local/projects/EBMAL/


	for j in $( \ls /local/projects/EBMAL/ | grep "Bp1AM" | grep -v "FR3_UWO");
	do
		fastqnumber=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | wc -l )
		if [[ $fastqnumber -eq 2 ]]
		then
			fastq1=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | head -1)
			fastq2=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | tail -1)

			qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /home/jmattick/SNP_Scripts/Pahangi_FASTQ_Combine_EBMAL.sh -fastq1 $fastq1 -fastq2 $fastq2 -Sample $j > $outdir/$j.Combine.txt

			list=$(cat $outdir/$j.Combine.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			output=$(echo $j"_R1.fastq.gz")

			fastq1=$(echo /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output)
			output=$(echo $j"_R2.fastq.gz")

			fastq2=$(echo /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output)



			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

			list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
			
			list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

			list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf


		else
			fastq1=$(\ls /local/projects/EBMAL/$j/ILLUMINA_DATA/*_R1.fastq.gz | grep "fastq")
			fastq2=$(echo $fastq1 | sed 's/R1.fastq.gz/R2.fastq.gz/g')


			qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

			list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
			
			list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true > $outdir/$j.BPDedup.txt

			list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

			qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf

		fi
	done

fi




####Recall Variants
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/MultiAF/")
	cd $outdir
	for j in $( \ls $outdir/*dedup*bam);
	do
		vcf=$(echo $j | sed 's/bam/vcf/g')
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $j --output $vcf
	done

	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/Clinical/")
	cd $outdir
	for j in $( \ls $outdir/*dedup*bam);
	do
		vcf=$(echo $j | sed 's/bam/vcf/g')
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $j --output $vcf
	done

	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/FR3_old/")
	cd $outdir
	for j in $( \ls $outdir/*dedup*bam);
	do
		vcf=$(echo $j | sed 's/bam/vcf/g')
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $j --output $vcf
	done

	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/FR3_UWO/")
	cd $outdir
	for j in $( \ls $outdir/*dedup*bam);
	do
		vcf=$(echo $j | sed 's/bam/vcf/g')
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $j --output $vcf
	done

	outdir=$(echo "/local/scratch/jmattick/PahangiIllumina/Bp1AM/")
	cd $outdir
	for j in $( \ls $outdir/*dedup*bam);
	do
		vcf=$(echo $j | sed 's/bam/vcf/g')
		qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $j --output $vcf
	done

fi


if [[ 1 -eq 2 ]]
then
	rm -f /local/scratch/jmattick/PahangiIllumina/MetricsTable.tsv
	cat /local/scratch/jmattick/PahangiIllumina/Bp1AM/Bp1AM_01.sorted.metrics | grep "READ_PAIR_OPTICAL_DUPLICATES" >> /local/scratch/jmattick/PahangiIllumina/MetricsTable.tsv

	for j in $( \ls /local/scratch/jmattick/PahangiIllumina/*/*.metrics | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
	do
		file=$(\ls /local/scratch/jmattick/PahangiIllumina/*/*.metrics | grep $j)
		cat $file | grep -v "#" | grep $j >> /local/scratch/jmattick/PahangiIllumina/MetricsTable.tsv

	done
fi


###This table is exported to R for analysis and plotting
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/PahangiIllumina/PahangiRTable/
	mkdir /local/scratch/jmattick/PahangiIllumina/PahangiRTable/

	for j in $( \ls /local/scratch/jmattick/PahangiIllumina/*/*.vcf | grep -v "Bp-.-FR3")
	do
		name=$(echo $j | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')

		paste -d '\t' <(cat $j | grep -v "#" | awk '{print $1"\t"$2}') <(cat $j | grep -v "#" | grep -o ".\/." | awk -F '/' '{print $1"\t"$2}') > /local/scratch/jmattick/PahangiIllumina/PahangiRTable/$name.tsv

	done

	cd /local/scratch/jmattick/PahangiIllumina/
	zip -r PahangiRTable.zip ./PahangiRTable/
	cp PahangiRTable.zip ~/
fi

if [[ 1 -eq 2 ]]
then
	cd /local/scratch/jmattick/PahangiIllumina/
	MalayiPath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa)

	/usr/local/packages/mummer/promer -p PahangivsMalayiCURRENT $MalayiPath $fastapath
	/usr/local/packages/mummer/show-coords -qlTHb PahangivsMalayiCURRENT.delta > PahangivsMalayiCURRENT.coords
fi

if [[ 1 -eq 2 ]]
then
	count=1
	for j in $( \ls /local/scratch/jmattick/PahangiIllumina/*/*.depth | grep -v "Bp-.-FR3")
	do
		name=$(echo $j | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		if [[ $count -eq 1 ]]
		then
			cat $j | awk -v name="$name" '{if(NR == 1){print "Contig\tPos\t"name"\n"$0} else {print $0}}' > /local/scratch/jmattick/PahangiIllumina/Combined.DepthTable.tsv
			count=2
		else
			paste -d '\t' <(cat /local/scratch/jmattick/PahangiIllumina/Combined.DepthTable.tsv) <(cat $j | awk -v name="$name" '{if(NR == 1){print name"\n"$3} else {print $3}}') > /local/scratch/jmattick/PahangiIllumina/temp.tsv
			mv /local/scratch/jmattick/PahangiIllumina/temp.tsv /local/scratch/jmattick/PahangiIllumina/Combined.DepthTable.tsv
		fi
	done
	gzip /local/scratch/jmattick/PahangiIllumina/Combined.DepthTable.tsv
fi

if [[ 1 -eq 2 ]]
then
	
	for j in $(\ls /local/scratch/jmattick/PahangiIllumina/*/*.vcf | grep "vcf" | grep -v "idx")
	do
		bgzip -c $j > $j.gz
		tabix -p vcf $j.gz
	done
fi


if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/PahangiSNPs/")
	rm -rf $outdir
	mkdir $outdir
	filelist=$(\ls /local/scratch/jmattick/PahangiIllumina/*/*.vcf.gz | grep "vcf.gz" | grep -v "idx" | tr '\n' ' ')
	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-merge -R "0/0" $filelist > $outdir/All.merged.ref.vcf


fi


##Filtering of all VCF files
if [[ 1 -eq 2 ]]
then
	rm -rf /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/
	mkdir /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/
	cp /local/scratch/jmattick/PahangiIllumina/*/*.vcf /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/
fi

if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/Filtered/")
	rm -rf $outdir
	mkdir $outdir
	for j in $(\ls /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/ | grep "vcf")
	do
		output=$(echo $j | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $j | sed 's/vcf/filteredonly.vcf/g')
		outputfilteredtsv=$(echo $j | sed 's/vcf/filteredonly.tsv/g')

		/usr/local/packages/gatk-4.0.4.0/gatk VariantFiltration -V /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/$j -R $fastapath -O $outdir"/"$output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"

		/usr/local/packages/gatk-4.0.4.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outdir"/"$outputfiltered -select 'vc.isNotFiltered()'

		paste -d '\t' <(cat $outdir"/"$outputfiltered | grep -v "#" | grep "Chr" | awk '{print $1"\t"$2}') <(cat $outdir"/"$outputfiltered | grep -v "#" | grep "Chr" | awk '{print $NF}' | sed 's/:.*//g' | awk -F '/' '{if($1 == $2) print "hom"; else print "het";}') > $outdir"/"$outputfilteredtsv

	done
fi


###BM VCF Files: /local/aberdeen2rw/julie/JM_dir/FilteredVCFs/
###BP VCF Files: /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/Filtered/


##Move to T3
if [[ 1 -eq 2 ]]
then
	rm -rf /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/
	mkdir /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/
	cd /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/
	/home/jdhotopp/bin/residues.pl /local/scratch/jmattick/BrugiaPahangi.FINAL.V4.fasta > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/Bp.v4.res.txt

	/home/jdhotopp/bin/residues.pl /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/Bm.res.txt

	for j in $(\ls /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/Filtered/*.tsv | grep "tsv")
	do
		cp $j /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/
	done

	for j in $(\ls /local/aberdeen2rw/julie/JM_dir/FilteredVCFs/*.tsv | grep "tsv")
	do
		cp $j /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/
	done
fi



###Tajima's D and relatedness calculations
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/MergedVCF/")
	rm -rf $outdir
	mkdir $outdir
	#cd /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/Filtered/
	#for j in $(\ls /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/Filtered/*.vcf | grep "filteredonly.vcf" | grep -v "idx")
	#do
	#	bgzip -c $j > $j.gz
	#	tabix -p vcf $j.gz
	#done


	cd $outdir
	filelist=$(\ls /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/Filtered/*.vcf.gz | grep "filteredonly.vcf.gz" | grep -v "\-FR3" | grep -v "idx" | tr '\n' ' ')
	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-merge -R "0/0" $filelist > $outdir/All.merged.ref.vcf

	filelistnoclin=$(\ls /local/aberdeen2rw/julie/JM_dir/BP_FilteredVCFs/Filtered/*.vcf.gz | grep "filteredonly.vcf.gz" | grep -v "\-FR3" | grep -v "clinical" | grep -v "idx" | tr '\n' ' ')
	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-merge -R "0/0" $filelistnoclin > $outdir/All.merged.ref.noclin.vcf


	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-stats $outdir/All.merged.ref.vcf --prefix All.merged.stats

	/usr/local/packages/plink/plink2 --pca --out BP.plink2 --vcf All.merged.ref.vcf --double-id --allow-extra-chr --freq
	/usr/local/packages/plink/plink2 --pca --out BP.plink2.noclin --vcf All.merged.ref.noclin.vcf --double-id --allow-extra-chr --freq


	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.ref.vcf --TajimaD 100000 --out All.merged.vcftools.tjd
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.ref.vcf --het --out All.merged.vcftools.het
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.ref.vcf --relatedness --out All.merged.vcftools.relatedness
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.ref.vcf --window-pi 10000 --out All.merged.vcftools.pi
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.ref.vcf --depth --out All.merged.vcftools.depth
fi

##Build Sequencing Stats Table
if [[ 1 -eq 2 ]]
then
	echo -e "Sample\tTotalReads\tMappedReads\tPercentMapped\tPercentDuplicates\tGenomeCoverage" > ~/BP.Mapping.Summary.tsv
	GenomeSize=97832585
	for j in $( \ls /local/scratch/jmattick/Pahangi_Output/*/*.bam | grep "bam" | grep "dedup");
	do
		Sample=$(echo $j | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		metrics=$(echo $j | sed 's/dedup.bam/metrics/g')
		flagstat=$(echo $j | sed 's/bam/flagstat/g')
		totalreads=$(cat $flagstat | head -1 | sed 's/ +.*//g')
		totalmapped=$(cat $flagstat | grep "mapped (" | sed 's/ +.*//g')
		percmapped=$(cat $flagstat | grep "mapped (" | sed 's/.*(//g' | sed 's/ .*//g')
		coverage=$(echo $totalmapped | awk -v GenomeSize="$GenomeSize" '{print 100*$1/GenomeSize}')

		echo -e $Sample"\t"$totalreads"\t"$totalmapped"\t"$percmapped"\t"$coverage >> ~/BP.Mapping.Summary.tsv
	done


fi






##Wuch Merge

if [[ 1 -eq 2 ]]
then
	
	for j in $(\ls /local/scratch/jmattick/Wuch_Output/*.vcf | grep "vcf" | grep -v "idx")
	do
		bgzip -c $j > $j.gz
		tabix -p vcf $j.gz
	done
fi


if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/Wuch_Output_Merge/")
	rm -rf $outdir
	mkdir $outdir
	filelist=$(\ls /local/scratch/jmattick/Wuch_Output/*.vcf.gz | grep "vcf.gz" | grep -v "idx" | tr '\n' ' ')
	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-merge -R "0/0" $filelist > $outdir/All.merged.ref.vcf


fi


##B. Timori SNPs
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/Timori_Mapping/")
	rm -rf $outdir
	mkdir $outdir
	#mv ~/brugia_timori.PRJEB4663.WBPS14.genomic.fa.gz /local/scratch/jmattick/
	#gzip -d /local/scratch/jmattick/brugia_timori.PRJEB4663.WBPS14.genomic.fa.gz
	fastapath=$(echo "/local/scratch/jmattick/brugia_timori.PRJEB4663.WBPS14.genomic.fa")
	/usr/local/packages/bwa/bin/bwa index $fastapath
	/usr/local/packages/samtools/bin/samtools faidx $fastapath
	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath


	fastq1=$(echo "/local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_1.fastq")
	fastq2=$(echo "/local/aberdeen2rw/julie/JM_dir/GenomeSizes/OtherNematodes/FASTQ/ERR346916_2.fastq")

	j="BTimori"
	qsub -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.bam -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/bwa/bin/bwa mem -M -a -t 32 -R "@RG'\t'ID:"$j"'\t'LB:"$j"'\t'SM:"$j $fastapath $fastq1 $fastq2 > $outdir/$j.BPMappingIds.txt

	list=$(cat $outdir/$j.BPMappingIds.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=10G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$outdir/$j.bam O=$outdir/$j.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPSortingIDs.txt
		
	list=$(cat $outdir/$j.BPSortingIDs.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=12G -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$outdir/$j.sorted.bam O=$outdir/$j.sorted.dedup.bam M=$outdir/$j.sorted.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true TMP_DIR=/local/scratch/jmattick/temp/ > $outdir/$j.BPDedup.txt

	list=$(cat $outdir/$j.BPDedup.txt | grep -o "job.*" | sed 's/job //g' | sed 's/ .*//g' | paste -s -d,)

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.flagstat -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools flagstat $outdir/$j.sorted.dedup.bam

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -o $outdir/$j.sorted.dedup.depth -l mem_free=4G -b y -cwd /usr/local/packages/samtools/bin/samtools  depth -aa -d 10000000000 $outdir/$j.sorted.dedup.bam

		#qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk AddOrReplaceReadGroups --INPUT $outdir/$j.sorted.dedup.bam --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM $j --OUTPUT $outdir/$j.sorted.dedup.readgroup.bam --CREATE_INDEX true

		##Works. Why does BWA fuck it up?
		##samtools view -H $outdir/$j.sorted.dedup.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' |  samtools reheader - $outdir/$j.sorted.dedup.bam > $outdir/$j.sorted.dedup.reheader.bam

	qsub -hold_jid $list -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/gatk-4.0.4.0/gatk HaplotypeCaller --read-filter MappingQualityReadFilter --reference $fastapath --input $outdir/$j.sorted.dedup.bam --output $outdir/$j.sorted.dedup.vcf

fi

if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/Timori_Mapping/")
	j="BTimori"
	cd $outdir

	output=$(echo $j.sorted.dedup.vcf | sed 's/vcf/filtered.vcf/g')
	outputfiltered=$(echo $j.sorted.dedup.vcf | sed 's/vcf/filteredonly.vcf/g')
	outputfilteredtsv=$(echo $j.sorted.dedup.vcf | sed 's/vcf/filteredonly.tsv/g')
	fastapath=$(echo "/local/scratch/jmattick/brugia_timori.PRJEB4663.WBPS14.genomic.fa")


	/usr/local/packages/gatk-4.0.4.0/gatk VariantFiltration -V $outdir/$j.sorted.dedup.vcf -R $fastapath -O $outdir"/"$output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"

	/usr/local/packages/gatk-4.0.4.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outdir"/"$outputfiltered -select 'vc.isNotFiltered()'


	paste -d '\t' <(cat $outdir"/"$outputfiltered | grep -v "#" | awk '{print $1"\t"$2}') <(cat $outdir"/"$outputfiltered | grep -v "#" | awk '{print $NF}' | sed 's/:.*//g' | awk -F '/' '{if($1 == $2) print "hom"; else print "het";}') > $outdir"/"$outputfilteredtsv

	MalayiPath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa)

	/usr/local/packages/mummer/promer -p BMvsBTCURRENT $MalayiPath $fastapath
	/usr/local/packages/mummer/show-coords -qlTHb BMvsBTCURRENT.delta > BMvsBTCURRENT.coords

fi


if [[ 1 -eq 2 ]]
then
	for j in $( \ls /local/scratch/jmattick/Pahangi_Output/*/*.vcf | grep "vcf" | grep "dedup");
	do
		statsout=$(echo $j | sed 's/vcf/stats/g')
		/usr/local/packages/bcftools/bin/bcftools stats $j > $statsout
	done
	for j in $( \ls /local/scratch/jmattick/Malayi_Output/*.vcf | grep "vcf" | grep "dedup");
	do
		statsout=$(echo $j | sed 's/vcf/stats/g')
		/usr/local/packages/bcftools/bin/bcftools stats $j > $statsout
	done

	for j in $( \ls /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/*.vcf | grep "vcf");
	do
		statsout=$(echo $j | sed 's/vcf/stats/g')
		/usr/local/packages/bcftools/bin/bcftools stats $j > $statsout
	done


fi

if [[ 1 -eq 2 ]]
then
	echo -e "Sample\tSNPs\tIndels\tTsTv" > ~/BP.Variant.Summary.tsv

	for j in $( \ls /local/scratch/jmattick/Pahangi_Output/*/*.stats | grep "stats" | grep "dedup");
	do
		Sample=$(echo $j | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		SNPNumber=$(cat $j | grep -v "#" | grep "SN" | grep -v "multi" | grep -E "SNP" | sed 's/.*://g' | awk '{print $1}')
		IndelNumber=$(cat $j | grep -v "#" | grep "SN" | grep -v "multi" | grep -E "indel" | sed 's/.*://g' | awk '{print $1}')
		tstv=$(cat $j | grep -v "#" | grep "TSTV" | awk '{print $5}')
		echo -e $Sample"\t"$SNPNumber"\t"$IndelNumber"\t"$tstv >> ~/BP.Variant.Summary.tsv
	done


	echo -e "Sample\tSNPs\tIndels\tTsTv" > ~/BM.Variant.Summary.tsv

	for j in $( \ls /local/scratch/jmattick/Malayi_Output/*.stats | grep "stats" | grep "dedup");
	do
		Sample=$(echo $j | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		SNPNumber=$(cat $j | grep -v "#" | grep "SN" | grep -v "multi" | grep -E "SNP" | sed 's/.*://g' | awk '{print $1}')
		IndelNumber=$(cat $j | grep -v "#" | grep "SN" | grep -v "multi" | grep -E "indel" | sed 's/.*://g' | awk '{print $1}')
		tstv=$(cat $j | grep -v "#" | grep "TSTV" | awk '{print $5}')
		echo -e $Sample"\t"$SNPNumber"\t"$IndelNumber"\t"$tstv >> ~/BM.Variant.Summary.tsv
	done

	#Thai Variants
	echo -e "Sample\tSNPs\tIndels\tTsTv" > ~/BM.Thai.Variant.Summary.tsv
	for j in $( \ls /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/*.stats | grep "stats" | grep "dedup");
	do
		Sample=$(echo $j | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
		SNPNumber=$(cat $j | grep -v "#" | grep "SN" | grep -v "multi" | grep -E "SNP" | sed 's/.*://g' | awk '{print $1}')
		IndelNumber=$(cat $j | grep -v "#" | grep "SN" | grep -v "multi" | grep -E "indel" | sed 's/.*://g' | awk '{print $1}')
		tstv=$(cat $j | grep -v "#" | grep "TSTV" | awk '{print $5}')
		echo -e $Sample"\t"$SNPNumber"\t"$IndelNumber"\t"$tstv >> ~/BM.Thai.Variant.Summary.tsv
	done


fi


###LD for B. malayi
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BrugiaMalayiSF2/
	mkdir /local/scratch/jmattick/BrugiaMalayiSF2/
	cd /local/scratch/jmattick/BrugiaMalayiSF2/
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/MergedVCF/")
	/usr/local/packages/bcftools/bin/bcftools query -f '%CHROM %POS[\t%DP\t%AD]\n' $outdir/All.merged.vcf > Bm.allsamples.bcfquery.freq
	cat Bm.allsamples.bcfquery.freq | grep "Chr" > Bm.allsamples.bcfquery.chrfreq

	/usr/local/packages/vcftools/bin/vcftools --vcf $outdir/All.merged.vcf --freq --out Bm.all.freq

	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/beagle/beagle.jar gt=$outdir/All.merged.ref.vcf out=All.merged.phased.vcf window=50000 overlap=3000 niterations=5
	zcat All.merged.phased.vcf.vcf.gz | grep "#" > All.merged.phased.filter.vcf
	zcat All.merged.phased.vcf.vcf.gz | grep -v "#" | awk '($10 != "0|0" || $11 != "0|0" || $12 != "0|0" || $13 != "0|0") && ($14 != "0|0" || $15 != "0|0" || $24 != "0|0" || $25 != "0|0") && ($16 != "0|0" || $17 != "0|0" || $18 != "0|0" || $19 != "0|0") && ($20 != "0|0" || $21 != "0|0" || $22 != "0|0" || $23 != "0|0") && ($26 != "0|0" || $27 != "0|0" || $28 != "0|0" || $29 != "0|0") && ($30 != "0|0" || $31 != "0|0" || $32 != "0|0" || $33 != "0|0" || $34 != "0|0" || $35 != "0|0") {print $0}' >> All.merged.phased.filter.vcf
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/vcftools/bin/vcftools --gzvcf All.merged.phased.vcf.vcf.gz  --hap-r2 --ld-window-bp 50000 --out ld_window_50000
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/vcftools/bin/vcftools --vcf All.merged.phased.filter.vcf  --hap-r2 --ld-window-bp 200000 --out ld_window_200kb_filtered

	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/vcftools/bin/vcftools --vcf All.merged.phased.filter.vcf  --hap-r2 --ld-window-bp 200000 --out ld_window_200kb_filtered
	echo -e "position\tx\tn\tfolded" > Bm.all.freq.chr.frq
	cat Bm.all.freq.frq | grep "ChrX" | awk '{print $2"\t"$3"\t"$4"\t0"}' >> Bm.all.freq.chr.frq
	export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64:$LD_LIBRARY_PATH
	/local/scratch/jmattick/miniconda3/bin/SweepFinder2 -sg 1000 Bm.all.freq.chr.frq Bm.all.freq.chr.SF

	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/vcftools/bin/vcftools --gzvcf All.merged.phased.vcf.vcf.gz --TajimaD 50000

	cat out.Tajima.D | grep "Chr" | grep -v "nan" > ~/BM.TajD.tsv


	####Create samples
	zcat All.merged.phased.vcf.vcf.gz | grep "#" | tail -1 | awk '{print $10"\n"$11"\n"$12"\n"$13}' > Lucknow.pops.txt
	zcat All.merged.phased.vcf.vcf.gz | grep "#" | tail -1 | awk '{print $14"\n"$15"\n"$24"\n"$25}' > TRS.pops.txt
	zcat All.merged.phased.vcf.vcf.gz | grep "#" | tail -1 | awk '{print $16"\n"$17"\n"$18"\n"$19}' > Thai.pops.txt
	zcat All.merged.phased.vcf.vcf.gz | grep "#" | tail -1 | awk '{print $20"\n"$21"\n"$22"\n"$23}' > FR3.pops.txt
	zcat All.merged.phased.vcf.vcf.gz | grep "#" | tail -1 | awk '{print $26"\n"$27"\n"$28"\n"$29}' > Liverpool.pops.txt
	zcat All.merged.phased.vcf.vcf.gz | grep "#" | tail -1 | awk '{print $30"\n"$31"\n"$32"\n"$33"\n"$34"\n"$35}' > WashU.pops.txt
	cat TRS.pops.txt FR3.pops.txt Liverpool.pops.txt WashU.pops.txt > Recent.FR3Lineage.txt
	cat Recent.FR3Lineage.txt Lucknow.pops.txt > FR3Lineage.txt

	/usr/local/packages/vcftools/bin/vcftools --gzvcf All.merged.phased.vcf.vcf.gz --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Lucknow.pops.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/TRS.pops.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Thai.pops.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/FR3.pops.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Liverpool.pops.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/WashU.pops.txt --fst-window-size 50000 --fst-window-step 50000 --out FST_BMalayi_All
	/usr/local/packages/vcftools/bin/vcftools --gzvcf All.merged.phased.vcf.vcf.gz --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Lucknow.pops.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Recent.FR3Lineage.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Thai.pops.txt --fst-window-size 50000 --fst-window-step 50000 --out FST_BMalayi_Lineage
	/usr/local/packages/vcftools/bin/vcftools --gzvcf All.merged.phased.vcf.vcf.gz --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/FR3Lineage.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Thai.pops.txt --fst-window-size 50000 --fst-window-step 50000 --out FST_BMalayi_Pair
	/usr/local/packages/vcftools/bin/vcftools --gzvcf All.merged.phased.vcf.vcf.gz --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/FR3Lineage.txt --weir-fst-pop /local/scratch/jmattick/BrugiaMalayiSF2/Thai.pops.txt --out FST_BMalayi_Pair_Ind

fi

if [[ 1 -eq 1 ]]
then
	cd /local/scratch/jmattick/BrugiaMalayiSF2/

	\ls /local/scratch/jmattick/BrugiaMalayiSF2/*.pops.txt > /local/scratch/jmattick/BrugiaMalayiSF2/Population.list
	start=1
	totallines=$(cat /local/scratch/jmattick/BrugiaMalayiSF2/Population.list | wc -l)
	while read line
	do
		name1=$(echo $line | awk -F '/' '{print $NF}' | sed 's/\..*//g')
		otherfile=$(cat /local/scratch/jmattick/BrugiaMalayiSF2/Population.list)
		if [[ $start -lt $totallines ]]
		then
			cat /local/scratch/jmattick/BrugiaMalayiSF2/Population.list | awk -v start="$start" 'NR > start {print $0}' > /local/scratch/jmattick/BrugiaMalayiSF2/Population.list2
			while read line2
			do
				name2=$(echo $line2 | awk -F '/' '{print $NF}' | sed 's/\..*//g')
				output=$(echo $name1"_"$name2)

				#/usr/local/packages/vcftools/bin/vcftools --gzvcf All.merged.phased.vcf.vcf.gz --weir-fst-pop $line --weir-fst-pop $line2 --fst-window-size 50000 --fst-window-step 50000 --out $output
				weightedFST=$(cat $output.log | grep "Weir and Cockerham weighted Fst estimate" | sed 's/Weir and Cockerham weighted Fst estimate: //g')
				echo -e $name1"\t"$name2"\t"$weightedFST >> FSTTable.tsv
			done < /local/scratch/jmattick/BrugiaMalayiSF2/Population.list2
		fi
		start=$(($start + 1))
	done < /local/scratch/jmattick/BrugiaMalayiSF2/Population.list

fi

###LD for B. pahangi
if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BrugiaPahangiSF2/
	mkdir /local/scratch/jmattick/BrugiaPahangiSF2/
	cd /local/scratch/jmattick/BrugiaPahangiSF2/

	cp /local/scratch/jmattick/BP_Kinship/All.BP.merged.beagle.vcf.gz /local/scratch/jmattick/BrugiaPahangiSF2/


	/local/scratch/jmattick/miniconda3/bin/SweepFinder2
	#zcat All.BP.merged.beagle.vcf.gz | grep "#" > All.BP.merged.beagle.filter.vcf
	#zcat All.BP.merged.beagle.vcf.gz | grep -v "#" | awk '($10 != "0|0" || $11 != "0|0" || $12 != "0|0" || $18 != "0|0" || $19 != "0|0" || $20 != "0|0" || $21 != "0|0" || $22 != "0|0") && ($13 != "0|0" || $14 != "0|0" || $15 != "0|0") {print $0}' >> All.BP.merged.beagle.filter.vcf
	zcat All.BP.merged.beagle.vcf.gz | grep "#" > All.BP.merged.beagle.hardfilter.vcf
	zcat All.BP.merged.beagle.vcf.gz | grep -v "#" | awk '($10 != "0|0" || $11 != "0|0" || $12 != "0|0") && ($18 != "0|0" || $19 != "0|0" || $20 != "0|0" || $21 != "0|0") && ($22 != "0|0") && ($13 != "0|0" || $14 != "0|0" || $15 != "0|0") {print $0}' >> All.BP.merged.beagle.hardfilter.vcf

	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/vcftools/bin/vcftools --vcf All.BP.merged.beagle.filter.vcf  --hap-r2 --ld-window-bp 200000 --out ld_window_BP_200kb_filtered
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/vcftools/bin/vcftools --vcf All.BP.merged.beagle.hardfilter.vcf  --hap-r2 --ld-window-bp 200000 --out ld_window_BP_200kb_hardfiltered

	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=24G -b y -cwd /usr/local/packages/vcftools/bin/vcftools --gzvcf All.BP.merged.beagle.vcf.gz --TajimaD 50000

	cat out.Tajima.D | grep "Chr" | grep -v "nan" > ~/BP.TajD.tsv

	##FST
	#~/BP.SamplePops.tsv
	cat ~/BP.SamplePops.tsv | awk '$2 == "FR3" {print $1}' > /local/scratch/jmattick/BrugiaPahangiSF2/Fr3.pops.txt
	cat ~/BP.SamplePops.tsv | awk '$2 == "Endemic" {print $1}' > /local/scratch/jmattick/BrugiaPahangiSF2/Endemic.pops.txt

	/usr/local/packages/vcftools/bin/vcftools --gzvcf All.BP.merged.beagle.vcf.gz --weir-fst-pop /local/scratch/jmattick/BrugiaPahangiSF2/Fr3.pops.txt --weir-fst-pop /local/scratch/jmattick/BrugiaPahangiSF2/Endemic.pops.txt --fst-window-size 50000 --fst-window-step 50000 --out FST_FR3_vs_Endemic

fi





