

###########################
# Input data
###########################

# List of non-redundent summary statistics
non_redundent_summary_statistics_file="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"

# Directory containing summary statistics
sumstat_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/sumstats/"

# Directory of quasi-independent LD blocks (from cTWAS package)
# Found here: https://github.com/xinhe-lab/ctwas/tree/main/inst/extdata/ldetect on April 3 2025
quasi_independent_ld_blocks_file="/lab-share/CHIP-Strober-e2/Public/LD_blocks/ldetect/EUR.b38.bed"

# Directory containing plink files for 100G
kg_plink_dir="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/"

# Directory containing baselineLD annotation
baselineLD_anno_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

# Directory containing ldsc weights
ldsc_weights_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

# File containing hapmap3 rsid
hapmap3_rsid_file="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"

# GTEx gencode gene annotation file
# wget https://storage.googleapis.com/adult-gtex/references/v8/reference-tables/gencode.v26.GRCh38.genes.gtf on 8/25/26
gene_annotation_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v26.GRCh38.genes.gtf"

# List of constrained genes
constrained_genes_geneset_file="/lab-share/CHIP-Strober-e2/Public/ben/snp_gene_disease_links/input_data/AMM_gs_constrained.txt"

# Directory containing pops results
# Downloaded from https://www.finucanelab.org/data on 6/20/23
pops_results_summary_file="/lab-share/CHIP-Strober-e2/Public/POPs/PoPS_FULLRESULTS.txt.gz"

# File containing MAGMA z-scores across traits
magma_z_score_file="/lab-share/CHIP-Strober-e2/Public/ben/snp_gene_disease_links/input_data/MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt"

# LDL silver standard gene sets
ldl_silver_standard_gene_set_file="/lab-share/CHIP-Strober-e2/Public/ben/snp_gene_disease_links/input_data/silver_standard_ldl_cholesterol_genes.csv"

# Names of GWAS traits
gwas_traits_file="/lab-share/CHIP-Strober-e2/Public/ben/snp_gene_disease_links/input_data/new_gwas_trait_names.txt"

# GWAS fine-mapping results
variant_fine_mapping_dir="/lab-share/CHIP-Strober-e2/Public/ben/snp_gene_disease_links/input_data/KUSHAL_FINEMAP_SCRIPT/RESULTS/"


###########################
# Output data
###########################
# Output root
output_root="/lab-share/CHIP-Strober-e2/Public/ben/snp_gene_disease_links/"

# Directory containing processed LD data
processed_ld_data_dir=${output_root}"processed_LD/"

# Directory containing snp-gene annotations
snp_gene_annotation_dir=${output_root}"snp_gene_annotations/"

# Pre-organized snp-gene annotations
preorganized_snp_gene_annotation_dir=${output_root}"preorganized_snp_gene_annotations/"

# Disease-specific tmp data dir
disease_specific_tmp_data_dir=${output_root}"disease_specific_tmp_data/"

# Learned snp-gene-disease links
learned_snp_gene_links_dir=${output_root}"learned_snp_gene_links/"

# Directory containing gene-set enrichment results
gene_set_enrichment_results_dir=${output_root}"gene_set_enrichment_analyses/"

# Directorying containing additive snp-gene links
additive_snp_gene_links_dir=${output_root}"additive_snp_gene_links/"

visualize_results_dir=${output_root}"visualize_results_dir/"


##############################
# Code
##############################
if false; then
sh process_LD_data.sh $hapmap3_rsid_file $baselineLD_anno_dir $kg_plink_dir $quasi_independent_ld_blocks_file $processed_ld_data_dir $gwas_traits_file $sumstat_dir
fi




K_closest_genes="10"
if false; then
sbatch generate_snp_gene_annotations.sh $baselineLD_anno_dir $gene_annotation_file $K_closest_genes $snp_gene_annotation_dir $processed_ld_data_dir $hapmap3_rsid_file $preorganized_snp_gene_annotation_dir
fi




## Run inference
if false; then
prior_choice="inverse_gamma_1e-16"
method_version="snp_gene_component_fixed_to_smart_init"
sed 1d $gwas_traits_file | while read trait_name; do



	trait_file=${sumstat_dir}${trait_name}".sumstats"
	input_window_summary_file=${preorganized_snp_gene_annotation_dir}"window_LD_summary_with_snp_gene_anno_v2.txt"
	gene_summary_file=${preorganized_snp_gene_annotation_dir}"gene_name_to_integer_mapping.txt"
	output_stem=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}
	sbatch run_snp_gene_disease_linking.sh $trait_name $trait_file $input_window_summary_file ${gene_summary_file} ${disease_specific_tmp_data_dir}${trait_name} ${output_stem} $prior_choice $method_version

done
fi


# Run fine-mapped additivie snp-gene linking
if false; then
sed 1d $gwas_traits_file | while read trait_name; do
	echo $trait_name
	additive_sgl_output_file=${additive_snp_gene_links_dir}${trait_name}"addative_snp_gene_links_gene_summary.txt"
	sbatch run_additive_snp_gene_linking.sh $trait_name ${gene_annotation_file} $variant_fine_mapping_dir $additive_sgl_output_file $baselineLD_anno_dir
done
fi


if false; then
prior_choice="inverse_gamma_1e-16"
method_version="snp_gene_component_fixed_to_smart_init"
#output_stem=${gene_set_enrichment_results_dir}${trait_name}"_"$prior_choice"_"${method_version}
#sgdlinks_gene_summary_file=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}"_lmm_snp_gene_link_"${prior_choice}"_"${method_version}"_gene_score_averaged.txt"
method_identifier=$prior_choice"_"${method_version}
sbatch gene_set_enrichment_analyses.sh $gwas_traits_file $learned_snp_gene_links_dir $pops_results_summary_file $magma_z_score_file $ldl_silver_standard_gene_set_file $gene_set_enrichment_results_dir $method_identifier $additive_snp_gene_links_dir
fi




################
# OLD
#################

prior_choice="inverse_gamma_cross_gene_prior_1e-1"
method_version="snp_gene_component"
prior_choice_arr=( "inverse_gamma_1e-16" "inverse_gamma_1e-8" "inverse_gamma_1e-3" "inverse_gamma_cross_gene_prior_1e-2" "inverse_gamma_cross_gene_prior_1e-1" "inverse_gamma_cross_gene_prior_1e-0" "inverse_gamma_cross_gene_prior_1e1" )


method_version="snp_gene_component_fixed_to_smart_init"
prior_choice_arr=( "inverse_gamma_1e-16"  )
if false; then
sed 1d $gwas_traits_file | while read trait_name; do

for prior_choice in "${prior_choice_arr[@]}"; do
echo $trait_name"_"${prior_choice}

trait_file=${sumstat_dir}${trait_name}".sumstats"
input_window_summary_file=${preorganized_snp_gene_annotation_dir}"window_LD_summary_with_snp_gene_anno_v2.txt"
gene_summary_file=${preorganized_snp_gene_annotation_dir}"gene_name_to_integer_mapping.txt"
output_stem=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}
sbatch run_snp_gene_disease_linking.sh $trait_name $trait_file $input_window_summary_file ${gene_summary_file} ${disease_specific_tmp_data_dir}${trait_name} ${output_stem} $prior_choice $method_version $constrained_genes_geneset_file

done
done
fi




prior_choice="inverse_gamma_1e-16"
#output_stem=${gene_set_enrichment_results_dir}${trait_name}"_"$prior_choice"_"${method_version}
#sgdlinks_gene_summary_file=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}"_lmm_snp_gene_link_"${prior_choice}"_"${method_version}"_gene_score_averaged.txt"
method_identifier=$prior_choice"_"${method_version}
if false; then
sh gene_set_enrichment_analyses.sh $gwas_traits_file $learned_snp_gene_links_dir $pops_results_summary_file $magma_z_score_file $ldl_silver_standard_gene_set_file $gene_set_enrichment_results_dir $method_identifier $additive_snp_gene_links_dir
fi



if false; then
module load R/3.5.1
Rscript visualize_snp_gene_disease_linking.R ${learned_snp_gene_links_dir} ${gene_set_enrichment_results_dir} ${visualize_results_dir}
fi

















trait_name="UKB_460K.body_HEIGHTz"
trait_name="UKB_460K.disease_AID_ALL"
trait_name="UKB_460K.biochemistry_LDLdirect"
trait_name="UKB_460K.body_HEIGHTz"
trait_name="UKB_460K.bp_SYSTOLICadjMEDz"
trait_name="UKB_460K.biochemistry_LDLdirect"


if false; then
method_version="snp_gene_component_fixed_to_smart_init"
prior_choice="inverse_gamma_cross_gene_prior_1e-1"
output_stem=${gene_set_enrichment_results_dir}${trait_name}"_"$prior_choice"_"${method_version}
sgdlinks_gene_summary_file=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}"_lmm_snp_gene_link_"${prior_choice}"_"${method_version}"_gene_score_averaged.txt"
sh gene_set_enrichment_analyses.sh $trait_name $sgdlinks_gene_summary_file $pops_results_summary_file $magma_z_score_file $ldl_silver_standard_gene_set_file $output_stem
fi


prior_choice="inverse_gamma_cross_gene_prior_1e-0"
method_version="snp_gene_component_fixed_to_smart_init"
output_stem=${gene_set_enrichment_results_dir}${trait_name}"_"$prior_choice"_"${method_version}
sgdlinks_gene_summary_file=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}"_lmm_snp_gene_link_"${prior_choice}"_"${method_version}"_gene_score_averaged.txt"
if false; then
sh gene_set_enrichment_analyses.sh $trait_name $sgdlinks_gene_summary_file $pops_results_summary_file $magma_z_score_file $ldl_silver_standard_gene_set_file $output_stem
fi

if false; then

prior_choice="inverse_gamma_cross_gene_prior_1e-1"
method_version="snp_gene_component_fixed_to_smart_init"
output_stem=${gene_set_enrichment_results_dir}${trait_name}"_"$prior_choice"_"${method_version}
sgdlinks_gene_summary_file=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}"_lmm_snp_gene_link_"${prior_choice}"_"${method_version}"_gene_score_averaged.txt"
echo $prior_choice
sh gene_set_enrichment_analyses.sh $trait_name $sgdlinks_gene_summary_file $pops_results_summary_file $magma_z_score_file $ldl_silver_standard_gene_set_file $output_stem


prior_choice="inverse_gamma_cross_gene_prior_1e-2"
method_version="snp_gene_component_fixed_to_smart_init"
output_stem=${gene_set_enrichment_results_dir}${trait_name}"_"$prior_choice"_"${method_version}
sgdlinks_gene_summary_file=${learned_snp_gene_links_dir}"snp_gene_links_"${trait_name}"_lmm_snp_gene_link_"${prior_choice}"_"${method_version}"_gene_score_averaged.txt"
echo $prior_choice
sh gene_set_enrichment_analyses.sh $trait_name $sgdlinks_gene_summary_file $pops_results_summary_file $magma_z_score_file $ldl_silver_standard_gene_set_file $output_stem
fi



if false; then
module load R/3.5.1
fi
if false; then
Rscript visualize_snp_gene_disease_linking.R ${output_stem} ${visualize_results_dir}
fi









########
# OLD

if false; then
prior_choice="inverse_gamma_1e-3"
prior_choice="inverse_gamma_1e-16"
prior_choice="inverse_gamma_1e-8"
prior_choice="inverse_gamma_cross_gene_prior_1e-4"

method_version="null_component"
method_version="snp_gene_component"



prior_choice_arr=( "inverse_gamma_1e-16" "inverse_gamma_1e-10" "inverse_gamma_1e-7" "inverse_gamma_1e-3" "inverse_gamma_cross_gene_prior_1e-4" "inverse_gamma_cross_gene_prior_1e-3" "inverse_gamma_cross_gene_prior_1e-2" "inverse_gamma_cross_gene_prior_1e-1" "inverse_gamma_cross_gene_prior_1e-0")
method_version_arr=( "null_component" "snp_gene_component")
################################

# Loop through methods
for prior_choice in "${prior_choice_arr[@]}"; do
for method_version in "${method_version_arr[@]}"; do

sbatch run_snp_gene_disease_linking.sh $trait_name $trait_file $input_window_summary_file ${gene_summary_file} ${disease_specific_tmp_data_dir}${trait_name} ${output_stem} $prior_choice $method_version


done
done
fi




