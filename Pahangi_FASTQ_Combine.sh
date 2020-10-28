for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-SRA1" ]] 
		then
		var=$(($var+1))
		SRA1=${!var}
		echo "Variable SRA1 is "$SRA1
	elif [[ "${!var}" = "-SRA2" ]] 
		then
		var=$(($var+1))
		SRA2=${!var}
		echo "Variable SRA2 is "$SRA2
	elif [[ "${!var}" = "-Sample" ]] 
		then
		var=$(($var+1))
		Sample=${!var}
		echo "Variable Sample is "$Sample
	else
		echo "Unrecognized input, exitting"
		exit 1
	fi
done


if [[ $SRA2 = "none" ]] 
then
	output=$(echo $Sample"_1.fastq.gz")
	fastq1=$(echo $SRA1"_1.fastq.gz")
	cp /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_FASTQ/$fastq1 /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output

	output=$(echo $Sample"_2.fastq.gz")
	fastq1=$(echo $SRA1"_2.fastq.gz")
	cp /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_FASTQ/$fastq1 /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output

else
	output=$(echo $Sample"_1.fastq.gz")

	fastq1=$(echo $SRA1"_1.fastq.gz")
	fastq2=$(echo $SRA2"_1.fastq.gz")

	cat /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_FASTQ/$fastq1 /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_FASTQ/$fastq2 > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output

	output=$(echo $Sample"_2.fastq.gz")

	fastq1=$(echo $SRA1"_2.fastq.gz")
	fastq2=$(echo $SRA2"_2.fastq.gz")

	cat /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_FASTQ/$fastq1 /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_FASTQ/$fastq2 > /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_Pahangi_Final_FASTQ/$output

fi
