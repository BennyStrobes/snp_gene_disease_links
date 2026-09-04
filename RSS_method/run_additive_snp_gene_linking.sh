#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=10G                        # Memory total in MiB (for all cores)



trait_name="${1}"
gene_annotation_file="${2}"
variant_fine_mapping_dir="${3}"
additive_sgl_output_file="${4}"
baselineLD_anno_dir="${5}"




source ~/.bash_profile


python run_additive_snp_gene_linking.py $trait_name ${gene_annotation_file} $variant_fine_mapping_dir $additive_sgl_output_file $baselineLD_anno_dir