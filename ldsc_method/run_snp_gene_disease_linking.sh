#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)



trait_name="${1}"
trait_sumstat_file="${2}"
input_window_summary_file="${3}"
gene_summary_file="${4}"
output_stem="${5}"


conda deactivate
module purge
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_gpu/bin/activate


python3 run_snp_gene_disease_linking.py $trait_name $trait_sumstat_file $input_window_summary_file $gene_summary_file ${output_stem}