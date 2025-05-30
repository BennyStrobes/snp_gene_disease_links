#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)




hapmap3_rsid_file="${1}"
baselineLD_anno_dir="${2}"
kg_plink_dir="${3}"
quasi_independent_ld_blocks_file="${4}"
processed_ld_data_dir="${5}"


genome_wide_ld_windows_file=${processed_ld_data_dir}"genome_wide_ld_windows.txt"
python3 generate_genome_wide_ld_windows_file.py $baselineLD_anno_dir $quasi_independent_ld_blocks_file $genome_wide_ld_windows_file




snp_set="hampmap3_snps"
module load python/3.7.4
for chrom_num in $(seq 1 22); do 
	python3 generate_ld_matrices_for_ld_windows_on_single_chromosome.py $chrom_num $genome_wide_ld_windows_file $hapmap3_rsid_file $baselineLD_anno_dir $kg_plink_dir $snp_set $processed_ld_data_dir
done





python3 merge_ld_matrix_files_across_chromosomes.py  $processed_ld_data_dir $snp_set
