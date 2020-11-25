#!/bin/bash

#file for automated search of guide+pam in reference and variant genomes

ref_folder=$(realpath $1)
vcf_folder=$(realpath $2)
guide_file=$(realpath $3)
pam_file=$(realpath $4)
annotation_file=$(realpath $5)
sampleID=$(realpath $6)

bMax=$7
mm=$8
bDNA=$9
bRNA=${10}

merge_t=${11}

output_folder=$(realpath ${12})

touch $output_folder/log.txt

echo "Complete search Start:" $(date +%F-%T) >> $output_folder/log.txt

starting_dir=${13}
ncpus=${14}

echo "CPU used: $ncpus" 
#echo $ref_folder
#echo $vcf_folder
#echo $guide_file
#echo $pam_file
#echo $annotation_file
#echo $output_folder
ref_name=$(basename $1)
#folder_of_folders=$(dirname $1)
vcf_name=$(basename $2)
guide_name=$(basename $3)
pam_name=$(basename $4)
annotation_name=$(basename $5)

declare -a real_chroms
for file_chr in "$ref_folder"/*.fa
do
	file_name=$(basename $file_chr)
	chr=$(echo $file_name | cut -f 1 -d'.')
	echo "$chr"
	real_chroms+=("$chr")	
done

if [ "$vcf_name" != "_" ]; then
	declare -a array_fake_chroms
	for file_chr in "$vcf_folder"/*.vcf.gz
	do
		file_name=$(basename $file_chr)
		chr=$(echo $file_name | cut -f 2 -d'.')
		echo "fake$chr"
		array_fake_chroms+=("fake$chr")	
	done
fi

if ! [ -d "$output_folder" ]; then
	mkdir "$output_folder"
fi
cd "$output_folder/"

fullseqpam=$(cut -f1 -d' ' "$pam_file")
pos=$(cut -f2 -d' ' "$pam_file")
if [ $pos -gt 0 ]; then
	true_pam=${fullseqpam:${#fullseqpam}-$pos}
else
	true_pam=${fullseqpam:0:-$pos}
fi

cd "$starting_dir"

echo "Start post-analysis"

if [ "$vcf_name" != "_" ]; then
	dict_folder="$output_folder/dictionaries_$vcf_name/"
	echo "Post-analysis SNPs Start: "$(date +%F-%T) >> $output_folder/log.txt	
	final_res="$output_folder/final_results_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
	final_res_alt="$output_folder/final_results_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
	if ! [ -f "$final_res" ]; then
		touch "$final_res"
	else
		echo "This result already exists. $final_res"
		exit
	fi
	if ! [ -f "$final_res_alt" ]; then
		touch "$final_res_alt"
	else
		echo "This result already exists. $final_res_alt"
		exit
	fi
	
	./pool_post_analisi_snp.py $output_folder $ref_folder $vcf_name $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $sampleID $dict_folder $final_res $final_res_alt
	
	echo "Post-analysis SNPs End: "$(date +%F-%T) >> $output_folder/log.txt	

else

	final_res="$output_folder/final_results_${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
	final_res_alt="$output_folder/final_results_${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
	if ! [ -f "$final_res" ]; then
		touch "$final_res"
	else
		echo "This result already exists. $final_res"
		exit
	fi
	if ! [ -f "$final_res_alt" ]; then
		touch "$final_res_alt"
	else
		echo "This result already exists. $final_res_alt"
		exit
	fi
	for key in "${!real_chroms[@]}"
	do
		echo "Processing $key"
		LC_ALL=C grep -P "$key\t" "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		touch "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		./scriptAnalisiNNN_v3.sh "$output_folder/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key" "$annotation_file" "_" "$ref_folder" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder"
		rm "$output_folder/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		header=$(head -1 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt")
		tail -n +2 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >> "$final_res" #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		tail -n +2 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >> "$final_res_alt" #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
		rm "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
		rm "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"
		
	done

fi

echo "Adding header to files"
sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref" "$final_res"
sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref" "$final_res_alt"

if [ "$vcf_name" != "_" ]; then
	echo "SNPs analysis ended. Starting INDELs analysis"
	cd "$starting_dir"
	
	echo "Post-analysis INDELs Start: "$(date +%F-%T) >> $output_folder/log.txt	
	./pool_post_analisi_indel.py $output_folder $ref_folder $vcf_folder $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $sampleID "$output_folder/log_indels_$vcf_name" $final_res $final_res_alt
	echo "Post-analysis INDELs End: "$(date +%F-%T) >> $output_folder/log.txt	
	# for key in "${!array_fake_chroms[@]}"
	# do
		# echo "Checking for INDELs results for $key"
		# if ! [ -f "$output_folder/finished$key.txt" ]; then 
			# continue
		# else
			# echo "Found results for $key, starting post-analysis"
			# true_chr=$key
			# fake_chr="fake$true_chr"
			
			# LC_ALL=C fgrep $true_chr "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
			# ./analisi_indels_NNN.sh "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr" "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}" "$annotation_file" "$output_folder/log_indels_$vcf_name" "$ref_folder/$true_chr.fa" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder" 
			# rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
			# tail -n +2 "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt" >> "$final_res" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
			# tail -n +2 "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt" >> "$final_res_alt" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
			# rm "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
			# rm "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
			# rm "$output_folder/finished$key.txt"
		# fi
	# done

fi



echo "Adjusting Results Start: "$(date +%F-%T) >> $output_folder/log.txt	
cd "$starting_dir"
echo "Adjusting final results - sorting"
header=$(head -1 $final_res)
tail -n +2 "$final_res" | LC_ALL=C sort -k5,5 -k7,7n -o "$final_res.sorted"
sed -i 1i"$header" "$final_res.sorted"
mv "$final_res.sorted" "$final_res"
if [ "$vcf_name" != "_" ]; then
	header=$(head -1 $final_res_alt)
	tail -n +2 "$final_res_alt" | LC_ALL=C sort -k5,5 -k7,7n -o "$final_res_alt.sorted"
	sed -i 1i"$header" "$final_res_alt.sorted"
	mv "$final_res_alt.sorted" "$final_res_alt"
fi
echo "Adjusting Results End: "$(date +%F-%T) >> $output_folder/log.txt	

echo "Merging Close Targets Start: "$(date +%F-%T) >> $output_folder/log.txt	
./merge_close_targets_cfd.sh $final_res $final_res.trimmed $merge_t
rm $final_res

./merge_close_targets_cfd.sh $final_res_alt $final_res_alt.trimmed $merge_t
rm $final_res_alt
echo "Merging Close Targets End: "$(date +%F-%T) >> $output_folder/log.txt	

echo "Merging Alternative Chromosomes Start: "$(date +%F-%T) >> $output_folder/log.txt	
./merge_alt_chr.sh $final_res.trimmed $final_res.trimmed.chr_merged
rm $final_res.trimmed

./merge_alt_chr.sh $final_res_alt.trimmed $final_res_alt.trimmed.chr_merged
rm $final_res_alt.trimmed
echo "Merging Alternative Chromosomes End: "$(date +%F-%T) >> $output_folder/log.txt	

mv $final_res.trimmed.chr_merged $final_res
mv $final_res_alt.trimmed.chr_merged $final_res_alt

sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res_alt"


echo "Post-analysis End: "$(date +%F-%T) >> $output_folder/log.txt
echo "POST-ANALYSIS END"

