#!/bin/bash

output_folder=$1
ref_folder=$2
ref_name=$(basename $2)
vcf_name=$3
guide_file=$4
guide_name=$(basename $4)
mm=$5
bDNA=$6
bRNA=$7
annotation_file=$8
annotation_name=$(basename $8)
pam_file=$9
pam_name=$(basename $9)
sampleID=${10}
dict_folder=${11}

final_res=${12}
final_res_alt=${13}

key=${14}

echo "Processing SNPs for $key"
LC_ALL=C grep -P "$key\t" "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
LC_ALL=C grep -P "$key\t" "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
./scriptAnalisiNNN_v3.sh "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key" "$annotation_file" "$dict_folder" "$ref_folder" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder"
rm "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
rm "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
# header=$(head -1 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt")
# tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >> "$final_res" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
# tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >> "$final_res_alt" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
# rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
# rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"