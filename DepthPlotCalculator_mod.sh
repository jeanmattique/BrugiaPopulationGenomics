for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-depth" ]] 
		then
		var=$(($var+1))
		depth=${!var}
		echo "Variable depth is "$depth
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

species=$(echo $depth | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
sample=$(echo $depth | awk -F '/' '{print $NF}' | awk -F '.' '{print $2}')

Nigonfile=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/$species/$species.PIGraphing.tsv | grep "tsv")
resfile=$(\ls /local/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/$species/$species.res.txt | grep "res.txt")


genomelength=$(cat $resfile | awk '{print $3}' | tail -1)
echo $genomelength
invalidlength=$(cat $resfile | grep -o "other:.*" | sed 's/other://g' | awk '{sum += $1} END {print sum}')
echo $invalidlength
truelength=$(cat $resfile | grep -o "other:.*" | sed 's/other://g' | awk -v genomelength="$genomelength" '{sum += $1} END {print genomelength-sum}')
echo $truelength


meandepth=$(cat $depth | awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')

depth_ten_percent=$(cat $depth | awk '$3 > 9 {print $0}' | wc -l | awk -v genomelength="$genomelength" '{print $0/genomelength}')
depth_ten_percent_mod=$(cat $depth | awk '$3 > 9 {print $0}' | wc -l | awk -v truelength="$truelength" '{print $0/truelength}')

coverage=$(cat $depth | awk '$3 > 0 {print $0}' | wc -l | awk -v genomelength="$genomelength" '{print $0/genomelength}')
coverage_mod=$(cat $depth | awk '$3 > 0 {print $0}' | wc -l | awk -v truelength="$truelength" '{print $0/truelength}')

echo -e $species"\t"$sample"\t"$genomelength"\t"$invalidlength"\t"$meandepth"\t"$coverage"\t"$coverage_mod"\t"$depth_ten_percent"\t"$depth_ten_percent_mod > $outdir/$species.$sample.depthinfo.tsv