

###########################
# Input data
###########################

# List of non-redundent summary statistics
non_redundent_summary_statistics_file="/n/groups/price/ldsc/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"

# Directory containing summary statistics
sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2024/sumstats/"

# Directory containing plink files for 100G
kg_plink_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# Directory containing baselineLD annotation
baselineLD_anno_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

# Directory containing ldsc weights
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

# File containing hapmap3 rsid
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"



###########################
# Output data
###########################
# Output root
tmp_output_root="/n/scratch/users/b/bes710/snp_gene_disease_links/"

# Directory containing processed LD data
processed_ld_data_dir=${tmp_output_root}"processed_LD/"

# Directory containing snp-gene annotations
snp_gene_annotation_dir=${tmp_output_root}"snp_gene_annotations/"

# Pre-organized snp-gene annotations
preorganized_snp_gene_annotation_dir=${tmp_output_root}"preorganized_snp_gene_annotations/"

# Learned snp-gene-disease links
learned_snp_gene_links_dir=${tmp_output_root}"learned_snp_gene_links/"




##############################
# Code
##############################
if false; then
sh process_LD_data.sh $hapmap3_rsid_file $baselineLD_anno_dir $kg_plink_dir $processed_ld_data_dir
fi


K_closest_genes="10"
if false; then
sh generate_snp_gene_annotations.sh $baselineLD_anno_dir $gene_annotation_file $K_closest_genes $snp_gene_annotation_dir $processed_ld_data_dir $preorganized_snp_gene_annotation_dir
fi






trait_name="UKB_460K.biochemistry_LDLdirect"
trait_file=${sumstat_dir}${trait_name}".sumstats"
input_window_summary_file=${preorganized_snp_gene_annotation_dir}"window_LD_summary_with_snp_gene_anno_v2.txt"
gene_summary_file=${preorganized_snp_gene_annotation_dir}"gene_name_to_integer_mapping.txt"
output_stem=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}
if false; then
sbatch run_snp_gene_disease_linking.sh $trait_name $trait_file $input_window_summary_file ${gene_summary_file} ${output_stem}
fi






if false; then
conda deactivate
module purge
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_gpu/bin/activate
python3 temp.py
fi


