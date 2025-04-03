#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)




baselineLD_anno_dir="${1}"
gene_annotation_file="${2}"
K_closest_genes="${3}"
snp_gene_annotation_dir="${4}"
processed_ld_data_dir="${5}"
preorganized_snp_gene_annotation_dir="${6}"



if false; then
for chrom_num in $(seq 1 22); do 
	echo $chrom_num
	python3 identify_k_closest_genes_to_each_snp.py $baselineLD_anno_dir $gene_annotation_file $K_closest_genes $chrom_num $snp_gene_annotation_dir
done
fi

echo "A"
python3 generate_snp_gene_annotations.py $snp_gene_annotation_dir

echo "B"
python3 preorganize_snp_gene_annotations.py $snp_gene_annotation_dir $processed_ld_data_dir $preorganized_snp_gene_annotation_dir

echo "C"
python3 convert_gene_names_to_integers.py $preorganized_snp_gene_annotation_dir

