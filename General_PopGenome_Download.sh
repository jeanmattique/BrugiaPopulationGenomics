
for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-sra" ]] 
		then
		var=$(($var+1))
		sra=${!var}
		echo "Variable sra is "$sra
	elif [[ "${!var}" = "-species" ]] 
		then
		var=$(($var+1))
		species=${!var}
		echo "Variable species is "$species
	elif [[ "${!var}" = "-sample" ]] 
		then
		var=$(($var+1))
		sample=${!var}
		echo "Variable sample is "$sample
	elif [[ "${!var}" = "-outdir" ]] 
		then
		var=$(($var+1))
		outdir=${!var}
		echo "Variable outdir is "$outdir
	else
		echo "Unrecognized input, exitting"
		exit 1
	fi
done

cd $outdir/$species/$sample/
#/usr/local/packages/sratoolkit-2.9.0/bin/fastq-dump --split-files $j

/usr/local/packages/sratoolkit-2.9.0/bin/fastq-dump --split-3 $sra
