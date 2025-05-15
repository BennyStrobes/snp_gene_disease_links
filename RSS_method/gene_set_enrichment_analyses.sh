#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)



gwas_traits_file="${1}"
learned_snp_gene_links_dir="${2}"
pops_results_summary_file="${3}"
magma_z_score_file="${4}"
ldl_silver_standard_gene_set_file="${5}"
gene_set_enrichment_results_dir="${6}"
method_identifier="${7}"






python3 run_gene_set_enrichment_analysis.py $gwas_traits_file $learned_snp_gene_links_dir $pops_results_summary_file $magma_z_score_file $ldl_silver_standard_gene_set_file $gene_set_enrichment_results_dir $method_identifier