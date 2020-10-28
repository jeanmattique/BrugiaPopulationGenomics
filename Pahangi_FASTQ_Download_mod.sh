for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-i" ]] 
		then
		var=$(($var+1))
		i=${!var}
		echo "Variable i is "$i
	elif [[ "${!var}" = "-j" ]] 
		then
		var=$(($var+1))
		j=${!var}
		echo "Variable j is "$j
	else
		echo "Unrecognized input, exitting"
		exit 1
	fi
done


cd /local/scratch/jmattick/PahangiFASTQ/$i
/usr/local/packages/sratoolkit-2.9.0/bin/fastq-dump --split-files $j

