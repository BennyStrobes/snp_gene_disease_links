args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

make_ldl_distance_prior_barchart <- function(prior_file) {

	df <- read.table(prior_file, header=TRUE, sep="\t")

	pp<-ggplot(data=df, aes(x=prior_name, y=prior_mean)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=prior_mean_lb, ymax=prior_mean_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="", y="Prior Probability", fill="")  +
  		theme(axis.text.x = element_text(angle = 90, hjust = 1))

  	return(pp)
}


make_pops_se_barplot <- function(pops_enrichment_summary_file, informal_trait_name) {
	df <- read.table(pops_enrichment_summary_file, header=TRUE, sep="\t")
	df$n_genes = factor(df$n_genes, levels=c(5,10,20,50,100,200,500))
	df$method = factor(df$method, levels=c("Magma", "sgdLinks"))
	pp<-ggplot(data=df, aes(x=n_genes, y=POPS_mean, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=POPS_lb, ymax=POPS_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="Top N genes", y="Average POPS score", fill="",title=informal_trait_name) 

  }



learned_snp_gene_links_dir <- args[1]
gene_set_enrichment_results_dir <- args[2]
vis_dir <- args[3]


method_version="inverse_gamma_1e-16_snp_gene_component_fixed_to_smart_init"

trait_name="UKB_460K.body_HEIGHTz"
informal_trait_name="Height"
pops_enrichment_summary_file <- paste0(gene_set_enrichment_results_dir, trait_name, "_", method_version, "_zscore_pops_enrichments_avg_pops_scores_per_threshold.txt")
pops_se_barplot1 <- make_pops_se_barplot(pops_enrichment_summary_file, informal_trait_name)

trait_name="UKB_460K.disease_AID_ALL"
informal_trait_name="Autoimmune disease"
pops_enrichment_summary_file <- paste0(gene_set_enrichment_results_dir, trait_name, "_", method_version, "_zscore_pops_enrichments_avg_pops_scores_per_threshold.txt")
pops_se_barplot2 <- make_pops_se_barplot(pops_enrichment_summary_file, informal_trait_name)

trait_name="UKB_460K.biochemistry_LDLdirect"
informal_trait_name="LDL Cholesterol"
pops_enrichment_summary_file <- paste0(gene_set_enrichment_results_dir, trait_name, "_", method_version, "_zscore_pops_enrichments_avg_pops_scores_per_threshold.txt")
pops_se_barplot3 <- make_pops_se_barplot(pops_enrichment_summary_file, informal_trait_name)


trait_name="UKB_460K.bp_SYSTOLICadjMEDz"
informal_trait_name="Systolic BP"
pops_enrichment_summary_file <- paste0(gene_set_enrichment_results_dir, trait_name, "_", method_version, "_zscore_pops_enrichments_avg_pops_scores_per_threshold.txt")
pops_se_barplot4 <- make_pops_se_barplot(pops_enrichment_summary_file, informal_trait_name)


trait_name="UKB_460K.lung_FEV1FVCzSMOKE"
informal_trait_name="FEV1:FVC"
pops_enrichment_summary_file <- paste0(gene_set_enrichment_results_dir, trait_name, "_", method_version, "_zscore_pops_enrichments_avg_pops_scores_per_threshold.txt")
pops_se_barplot5 <- make_pops_se_barplot(pops_enrichment_summary_file, informal_trait_name)

legender <- get_legend(pops_se_barplot1)
joint_pops <- plot_grid(pops_se_barplot1 + theme(legend.position="none"), pops_se_barplot2+ theme(legend.position="none"), pops_se_barplot3+ theme(legend.position="none"), pops_se_barplot4+ theme(legend.position="none"), pops_se_barplot5+ theme(legend.position="none"), legender, ncol=2)



output_file <- paste0(vis_dir, "POPS_enrichment_barplot_", method_version, ".pdf")
ggsave(joint_pops, file=output_file, width=7.2, height=6.5, units="in")




if (FALSE) {
LDL_prior_file <- paste0(results_output_stem, "_lmm_snp_gene_link_inverse_gamma_cross_gene_prior_1e-1_snp_gene_component_snp_gene_priors_averaged.txt")


pp <- make_ldl_distance_prior_barchart(LDL_prior_file)
output_file <- paste0(vis_dir, "LDL_distance_prior_barchart.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")

}