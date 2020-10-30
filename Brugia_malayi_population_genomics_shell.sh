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

	/usr/local/packages/bwa/bin/bwa index $fastapath
	rm $dictpath
	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath

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

	done
fi


###Filter Variants
if [[ 1 -eq 2 ]]
then
	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa)
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/FilteredVCFs/")
	rm -rf $outdir
	mkdir $outdir
	for j in $( \ls /local/aberdeen2rw/julie/JM_dir/BMalayiCenters/ | grep "vcf" | grep -v "idx");
	do
		inputvcf=$(echo "/local/aberdeen2rw/julie/JM_dir/BMalayiCenters/"$j)
		output=$(echo $j | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $j | sed 's/vcf/filteredonly.vcf/g')

		#QD < 5, QUAL < 30, DP < 14, MQ < 30, MQRankSum < -12.5, ReadPosRankSum < -8.0, FS > 60.0, ABHet < .30, ABHet> 0.70, ABHom < 0.90

		/usr/local/packages/gatk-4.0.4.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $outdir"/"$output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"

		/usr/local/packages/gatk-4.0.4.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outdir"/"$outputfiltered -select 'vc.isNotFiltered()'

	done

	for j in $( \ls /local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/ | grep "vcf" | grep -v "idx");
	do
		inputvcf=$(echo "/local/aberdeen2rw/julie/JM_dir/BMalayiMalaysia/"$j)
		output=$(echo $j | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $j | sed 's/vcf/filteredonly.vcf/g')

		#QD < 5, QUAL < 30, DP < 14, MQ < 30, MQRankSum < -12.5, ReadPosRankSum < -8.0, FS > 60.0, ABHet < .30, ABHet> 0.70, ABHom < 0.90

		/usr/local/packages/gatk-4.0.4.0/gatk VariantFiltration -V $inputvcf -R $fastapath -O $outdir"/"$output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"

		/usr/local/packages/gatk-4.0.4.0/gatk SelectVariants -V $outdir"/"$output -R $fastapath -O $outdir"/"$outputfiltered -select 'vc.isNotFiltered()'

	done
fi



###Tajima's D and Relatedness Statistics
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/MergedVCF/")
	rm -rf $outdir
	mkdir $outdir
	cd /local/aberdeen2rw/julie/JM_dir/FilteredVCFs/
	for j in $(\ls /local/aberdeen2rw/julie/JM_dir/FilteredVCFs/*.vcf | grep "filteredonly.vcf" | grep -v "idx")
	do
		bgzip -c $j > $j.gz
		tabix -p vcf $j.gz
	done

	for j in $(\ls /local/aberdeen2rw/julie/JM_dir/FilteredVCFs/*.vcf | grep "filteredonly.vcf" | grep -v "idx")
	do
		paste -d '\t' <(cat $j | grep -v "#" | grep "Chr" | awk '{print $1"\t"$2}') <(cat $j | grep -v "#" | grep "Chr" | awk '{print $NF}' | sed 's/:.*//g' | awk -F '/' '{if($1 == $2) print "hom"; else print "het";}') > $j.tsv
	done


	cd $outdir
	filelist=$(\ls /local/aberdeen2rw/julie/JM_dir/FilteredVCFs/*.vcf.gz | grep "filteredonly.vcf.gz" | grep -v "idx" | tr '\n' ' ')
	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-merge $filelist > $outdir/All.merged.vcf
	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-merge -R "0/0" $filelist > $outdir/All.merged.ref.vcf

	perl -I /usr/local/packages/vcftools/lib/site_perl/5.24.0/ /usr/local/packages/vcftools/bin/vcf-stats $outdir/All.merged.vcf --prefix All.merged.stats

	/usr/local/packages/plink/plink2 --pca --out BM.plink2 --vcf ./All.merged.ref.vcf --double-id --allow-extra-chr --freq

	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.vcf --TajimaD 100000 --out All.merged.vcftools.tjd
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.vcf --het --out All.merged.vcftools.het
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.vcf --relatedness --out All.merged.vcftools.relatedness
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.vcf --window-pi 10000 --out All.merged.vcftools.pi
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.vcf --depth --out All.merged.vcftools.depth
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.vcf --out All.merged.vcftools.fst --weir-fst-pop ./Sample_Malaysia.txt --weir-fst-pop ./Sample_FR3.txt --weir-fst-pop ./Sample_TRS.txt --weir-fst-pop ./Sample_Lucknow.txt --weir-fst-pop ./Sample_Liverpool.txt --weir-fst-pop ./Sample_WashU.txt
	/usr/local/packages/vcftools/bin/vcftools --vcf All.merged.vcf --out All.merged.vcftools.fst.malaysia.FR3 --weir-fst-pop ./Sample_Malaysia.txt --weir-fst-pop ./Sample_FR3.txt
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





##Filter VCF for Runs of Homozygosity
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/aberdeen2rw/julie/JM_dir/MergedVCF/")
	cd $outdir
	cp ~/RunsofHomozygosity.tsv ./
	cat /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf | grep "#" > /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.RoH.vcf

	while read line
	do
		chr=$(echo $line | awk '{print $1}')
		start=$(echo $line | awk '{print $2}')
		end=$(echo $line | awk '{print $3}')
		cat /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf | grep $chr | awk -v start="$start" -v end="$end" '$2 > start && $2 < end {print $0}' >> /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.RoH.vcf
	done < /local/aberdeen2rw/julie/JM_dir/MergedVCF/RunsofHomozygosity.tsv
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


####Kinship for BPahangi
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/BP_Kinship/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	#PHASE THE DATA
	qsub -P jhotopp-gcid-proj4b-filariasis -l mem_free=5G -q threaded.q -pe thread 32 -b y -cwd /usr/local/packages/jdk-8u151/bin/java -Xmx10g -jar /usr/local/packages/beagle/beagle.jar gt=/local/scratch/jmattick/PahangiSNPs/All.merged.ref.vcf nthreads=32 ibd=true out=$outdir/All.BP.merged.beagle niterations=20

	##Fix sample IDs
	#zcat /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.vcf.gz | sed 's/BM3/TRS-male_3-M1/g' | sed 's/BM4/TRS-male_4-M1/g' | sed 's/Bm_Malaysia/Bm-Malaysia/g' | sed 's/F3_male/F3-male/g' | sed 's/_male_M1/-male-M1/g'  | sed 's/_m1/-m1/g' | sed 's/_M1/-M1/g' | sed 's/TRS_male/TRS-male/g' | sed 's/W_male/W-male/g' | grep -E "Chr|#" > /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf
	#/usr/local/packages/plink/plink2 --vcf $outdir/All.BP.merged.beagle.vcf.gz --make-king --homozyg --allow-extra-chr --out All.BP.merged.phased.plink
	#/usr/local/packages/plink/plink2 --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf --export 'ped' --allow-extra-chr --out All.merged.phased.plink
	
	##Generate ped and map files!
	#Bp1AM_01        Bp1AM_02        Bp1AM_03        Bp-1-clinical   Bp-2-clinical   Bp-3-clinical   Bp-5-FR3        Bp-6-FR3        FR3_UWO_Bp1AM_01        FR3_UWO_Bp1AM_05        FR3_UWO_Bp1AM_06        FR3_UWO_Bp1AM_09        BEI_FR3_Bpahangi_multi_AF_a
	zcat $outdir/All.BP.merged.beagle.vcf.gz | sed 's/Bp1AM_/Bp1AM-/g' | sed 's/FR3_UWO_/FR3-UWO-/g' | sed 's/BEI_FR3_Bpahangi_multi_AF_a/BEI-FR3-Bpahangi-multi-AF-a/g' > $outdir/All.BP.merged.sample.beagle.vcf
	/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf $outdir/All.BP.merged.sample.beagle.vcf --recode --allow-extra-chr --out All.BP.merged.phased.plink.1.9


	#####GET ROH FILE FROM R -- put in home directory and bring back
	cp ~/BP_RunsofHomozygosity.tsv $outdir
	cat $outdir/All.BP.merged.sample.beagle.vcf | grep "#" > $outdir/All.BP.merged.sample.beagle.ROH.vcf

	while read line
	do
		chr=$(echo $line | awk '{print $1}')
		start=$(echo $line | awk '{print $2}')
		end=$(echo $line | awk '{print $3}')
		cat $outdir/All.BP.merged.sample.beagle.vcf | grep $chr | awk -v start="$start" -v end="$end" '$2 > start && $2 < end {print $0}' >> $outdir/All.BP.merged.sample.beagle.ROH.vcf
	done < $outdir/BP_RunsofHomozygosity.tsv


	####Do Kinship on RoH regions only
	/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf $outdir/All.BP.merged.sample.beagle.ROH.vcf --allow-extra-chr --het --out All.BP.merged.beagle.sample.RoH
	#/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.vcf --het --allow-extra-chr --out All.merged.phased.plink.1.9
	#cat /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.ref.vcf | sed 's/BM3/TRS-male_3-M1/g' | sed 's/BM4/TRS-male_4-M1/g' | sed 's/Bm_Malaysia/Bm-Malaysia/g' | sed 's/F3_male/F3-male/g' | sed 's/_male_M1/-male-M1/g'  | sed 's/_m1/-m1/g' | sed 's/_M1/-M1/g' | sed 's/TRS_male/TRS-male/g' | sed 's/W_male/W-male/g' | grep -E "Chr|#" > /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.ref.header.vcf

	#/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.ref.header.vcf --het --allow-extra-chr --out All.merged.chr

	#/usr/local/packages/plink-1.90.beta-3.6/bin/plink --vcf /local/aberdeen2rw/julie/JM_dir/MergedVCF/All.merged.beagle.sample.RoH.noX.vcf --allow-extra-chr --het --out All.merged.beagle.sample.RoH.noX
	cat $outdir/All.BP.merged.beagle.sample.RoH.het | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ~/All.BP.merged.beagle.sample.RoH.het
fi



if [[ 1 -eq 2 ]]
then
	rm -rf /local/scratch/jmattick/BMalayi_PacBio/
	mkdir /local/scratch/jmattick/BMalayi_PacBio/
	cd /local/scratch/jmattick/BMalayi_PacBio/
	for SRA in $(cat ~/B.malayi.PacBioSRA.txt | grep "SRR")
	do
		qsub -P jhotopp-gcid-proj4b-filariasis -N SRADownload_$SRA -l mem_free=12G -b y -cwd /home/jmattick/SNP_Scripts/BM_PacBio_download.sh -SRA $SRA
	done

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
	outdir=$(echo "/local/scratch/jmattick/BrugiaPopGenome_CombinedFastqs/")
	#rm -rf $outdir
	#mkdir $outdir
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
	outdir=$(echo "/local/scratch/jmattick/BrugiaPopGenome_CombinedFastqs/")

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
	outdir=$(echo "/local/scratch/jmattick/WuchPopGenome_Fastqs/")
	#rm -rf $outdir
	#mkdir $outdir

	for i in $(cat ~/Wuchereria.bancrofti.populations.SRA.txt | awk '{print $2}' | sort | uniq)
	do
		mkdir /local/scratch/jmattick/WuchPopGenome_Fastqs/$i/
		cd /local/scratch/jmattick/WuchPopGenome_Fastqs/$i/
		for j in $(cat ~/Wuchereria.bancrofti.populations.SRA.txt | grep $i | awk '{print $1}' | sort | uniq)
		do

			qsub -V -P jhotopp-gcid-proj4b-filariasis -l mem_free=4G -b y -cwd /home/jmattick/SNP_Scripts/WB_PopGenome_Download.sh -i $i -j $j

		done
		filename=$(echo $i | sed 's/\r//g')
		mv $outdir/$i $outdir/$filename

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
if [[ 1 -eq 1 ]]
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

	for fastq1 in $(\ls /local/scratch/jmattick/BrugiaPopGenome_CombinedFastqs/*_1.fastq.gz | grep "fastq");
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

	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS14.genomic_masked.fa)

	dictpath=$(echo $fastapath | sed 's/fa/dict/g')

	/usr/local/packages/bwa/bin/bwa index $fastapath
	samtools faidx $fastapath
	rm $dictpath
	/usr/local/packages/jdk-8u151/bin/java -Xmx12g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$fastapath O=$dictpath


	outdir=$(echo "/local/scratch/jmattick/Wuch_PopGenomeRemap/")

	rm -rf $outdir
	mkdir $outdir
	cd $outdir

	for fastq1 in $(\ls /local/scratch/jmattick/WuchPopGenome_CombinedFastqs/*_1.fastq.gz | grep "fastq");
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


###BIG FINAL FILTERING
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/local/scratch/jmattick/FilteredVariants/")
	rm -rf $outdir
	mkdir $outdir

	species="BMalayi"
	mkdir $outdir"/"$species

	fastapath=$(\ls /local/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa)



	for j in $(\ls /local/scratch/jmattick/Malayi_Output/ | grep "vcf" | grep -v "idx")
	do
		output=$(echo $j | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $j | sed 's/vcf/filteredonly.vcf/g')
		outputfilteredtsv=$(echo $j | sed 's/vcf/filteredonly.tsv/g')

		/usr/local/packages/gatk-4.0.4.0/gatk VariantFiltration -V /local/scratch/jmattick/Malayi_Output/$j -R $fastapath -O $outdir"/"$species"/"$output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"

		/usr/local/packages/gatk-4.0.4.0/gatk SelectVariants -V $outdir"/"$species"/"$output -R $fastapath -O $outdir"/"$species"/"$outputfiltered -select 'vc.isNotFiltered()'

		paste -d '\t' <(cat $outdir"/"$species"/"$outputfiltered | grep -v "#" | grep "Chr" | awk '{print $1"\t"$2}') <(cat $outdir"/"$species"/"$outputfiltered | grep -v "#" | grep "Chr" | awk '{print $NF}' | sed 's/:.*//g' | awk -F '/' '{if($1 == $2) print "hom"; else print "het";}') > $outdir"/"$species"/"$outputfilteredtsv

	done

	species="BPahangi"
	mkdir $outdir"/"$species

	fastapath=$(\ls /local/scratch/jmattick/BrugiaPahangi.FINAL.V5.2.masked.fasta)


	for j in $(\ls /local/scratch/jmattick/Pahangi_Output/*/*.vcf | grep "vcf" | grep -v "idx")
	do
		output=$(echo $j | awk -F '/' '{print $NF}' | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $j | awk -F '/' '{print $NF}' | sed 's/vcf/filteredonly.vcf/g')
		outputfilteredtsv=$(echo $j | awk -F '/' '{print $NF}' | sed 's/vcf/filteredonly.tsv/g')

		/usr/local/packages/gatk-4.0.4.0/gatk VariantFiltration -V $j -R $fastapath -O $outdir"/"$species"/"$output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"

		/usr/local/packages/gatk-4.0.4.0/gatk SelectVariants -V $outdir"/"$species"/"$output -R $fastapath -O $outdir"/"$species"/"$outputfiltered -select 'vc.isNotFiltered()'

		paste -d '\t' <(cat $outdir"/"$species"/"$outputfiltered | grep -v "#" | grep "Chr" | awk '{print $1"\t"$2}') <(cat $outdir"/"$species"/"$outputfiltered | grep -v "#" | grep "Chr" | awk '{print $NF}' | sed 's/:.*//g' | awk -F '/' '{if($1 == $2) print "hom"; else print "het";}') > $outdir"/"$species"/"$outputfilteredtsv

	done

	species="WBancrofti"
	mkdir $outdir"/"$species

	fastapath=$(\ls /local/scratch/jmattick/wuchereria_bancrofti.PRJNA275548.WBPS14.genomic_masked.fa)


	for j in $(\ls /local/scratch/jmattick/Wuch_Output/ | grep "vcf" | grep -v "idx")
	do
		output=$(echo $j | sed 's/vcf/filtered.vcf/g')
		outputfiltered=$(echo $j | sed 's/vcf/filteredonly.vcf/g')
		outputfilteredtsv=$(echo $j | sed 's/vcf/filteredonly.tsv/g')

		/usr/local/packages/gatk-4.0.4.0/gatk VariantFiltration -V /local/scratch/jmattick/Wuch_Output/$j -R $fastapath -O $outdir"/"$species"/"$output --filter-expression "QD < 5.0 || QUAL < 30.0 || DP < 14.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" --filter-name "HardFilter"

		/usr/local/packages/gatk-4.0.4.0/gatk SelectVariants -V $outdir"/"$species"/"$output -R $fastapath -O $outdir"/"$species"/"$outputfiltered -select 'vc.isNotFiltered()'

		paste -d '\t' <(cat $outdir"/"$species"/"$outputfiltered | grep -v "#" | awk '{print $1"\t"$2}') <(cat $outdir"/"$species"/"$outputfiltered | grep -v "#" | awk '{print $NF}' | sed 's/:.*//g' | awk -F '/' '{if($1 == $2) print "hom"; else print "het";}') > $outdir"/"$species"/"$outputfilteredtsv

	done

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

