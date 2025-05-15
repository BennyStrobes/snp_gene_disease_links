#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)



trait_name="${1}"
trait_sumstat_file="${2}"
input_window_summary_file="${3}"
gene_summary_file="${4}"
disease_specific_tmp_data_dir="${5}"
output_stem="${6}"
prior_choice="${7}"
method_version="${8}"
geneset_file="${9}"


if false; then
python3 get_disease_specific_data_ready_for_snp_gene_disease_linking.py $trait_sumstat_file $input_window_summary_file $gene_summary_file $disease_specific_tmp_data_dir
fi


new_input_window_summary_file=${disease_specific_tmp_data_dir}"_window_summary_file.txt"

echo ${trait_name}
echo ${prior_choice}"_"${method_version}
if false; then
python3 run_snp_gene_disease_linking.py $trait_name $new_input_window_summary_file $gene_summary_file ${output_stem} $prior_choice $method_version $geneset_file
fi

python3 organize_snp_gene_disease_link_results.py $trait_name ${output_stem} $prior_choice $method_version
