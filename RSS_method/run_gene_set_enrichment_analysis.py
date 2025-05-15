import numpy as np
import os
import sys
import pdb
import gzip


def create_dictionary_mapping_trait_names_to_pops_trait_names():
	dicti = {}
	dicti['UKB_460K.body_HEIGHTz'] = 'Height'
	dicti['UKB_460K.biochemistry_LDLdirect'] = 'LDLC'
	dicti['UKB_460K.lung_FEV1FVCzSMOKE'] = 'FEV1FVC'
	dicti['UKB_460K.disease_AID_ALL'] = 'AID_Combined'
	dicti['UKB_460K.bp_SYSTOLICadjMEDz'] = 'SBP'

	dicti['UKB_460K.biochemistry_Cholesterol'] = 'TC'
	dicti['UKB_460K.biochemistry_Glucose'] = 'Glucose'
	dicti['UKB_460K.biochemistry_HDLcholesterol'] = 'HDLC'
	dicti['UKB_460K.biochemistry_HbA1c'] = 'HbA1c'
	dicti['UKB_460K.biochemistry_Triglycerides'] = 'TG'
	dicti['UKB_460K.biochemistry_VitaminD'] = 'VitD'
	dicti['UKB_460K.blood_EOSINOPHIL_COUNT'] = 'Eosino'
	dicti['UKB_460K.blood_LYMPHOCYTE_COUNT'] = 'Lym'
	dicti['UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN'] = 'MCH'
	dicti['UKB_460K.blood_MONOCYTE_COUNT'] = 'Mono'
	dicti['UKB_460K.blood_PLATELET_COUNT'] = 'Plt'
	dicti['UKB_460K.blood_RED_COUNT'] = 'RBC'
	dicti['UKB_460K.blood_WHITE_COUNT'] = 'WBC'
	dicti['UKB_460K.bmd_HEEL_TSCOREz'] = 'eBMD'
	dicti['UKB_460K.body_BMIz'] = 'BMI'
	dicti['UKB_460K.body_WHRadjBMIz'] = 'WHRadjBMI'
	dicti['UKB_460K.bp_DIASTOLICadjMEDz'] = 'DBP'
	dicti['UKB_460K.disease_ASTHMA_DIAGNOSED'] = 'Asthma'
	dicti['UKB_460K.disease_HYPOTHYROIDISM_SELF_REP'] = 'Hypothyroidism'
	dicti['UKB_460K.disease_T2D'] = 'T2D'
	dicti['UKB_460K.other_MORNINGPERSON'] = 'Morning_Person'
	dicti['UKB_460K.repro_MENARCHE_AGE'] = 'Age_at_Menarche'

	'''
	dicti['CAD'] = 'CAD'
	dicti['TC'] = 'biochemistry_Cholesterol'
	dicti['Glucose'] = 'biochemistry_Glucose'
	dicti['HDLC'] = 'biochemistry_HDLcholesterol'
	dicti['HbA1c'] = 'biochemistry_HbA1c'
	dicti['LDLC'] = 'biochemistry_LDLdirect'
	dicti['TG'] = 'biochemistry_Triglycerides'
	dicti['VitD'] = 'biochemistry_VitaminD'
	dicti['Eosino'] = 'blood_EOSINOPHIL_COUNT'
	dicti['Lym'] = 'blood_LYMPHOCYTE_COUNT'
	dicti['MCH'] = 'blood_MEAN_CORPUSCULAR_HEMOGLOBIN'
	dicti['Mono'] = 'blood_MONOCYTE_COUNT'
	dicti['Plt'] = 'blood_PLATELET_COUNT'
	dicti['RBC'] = 'blood_RED_COUNT'
	dicti['WBC'] = 'blood_WHITE_COUNT'
	dicti['eBMD'] = 'bmd_HEEL_TSCOREz'
	dicti['BMI'] = 'body_BMIz'
	dicti['Height'] = 'body_HEIGHTz'
	dicti['WHRadjBMI'] = 'body_WHRadjBMIz'
	dicti['DBP'] = 'bp_DIASTOLICadjMEDz'
	dicti['AID_Combined'] = 'disease_AID_ALL'
	dicti['Asthma'] = 'disease_ASTHMA_DIAGNOSED'
	dicti['Hypothyroidism'] = 'disease_HYPOTHYROIDISM_SELF_REP'
	dicti['T2D'] = 'disease_T2D'
	dicti['Morning_Person'] = 'other_MORNINGPERSON'
	dicti['Age_at_Menarche'] = 'repro_MENARCHE_AGE'
	'''
	return dicti


def extract_pops_scores(pops_results_summary_file, pops_trait_name):
	ens_id_to_pops_score = {}
	gene_id_to_ens_id = {}
	f = gzip.open(pops_results_summary_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_name = data[0]
		if trait_name != pops_trait_name:
			continue
		ens_id = data[2]
		gene_id = data[3]
		pops_score = float(data[8])

		if ens_id in ens_id_to_pops_score:
			print('assumption eroro')
			pdb.set_trace()
		ens_id_to_pops_score[ens_id] = pops_score


		gene_id_to_ens_id[gene_id] = ens_id
	f.close()

	return ens_id_to_pops_score, gene_id_to_ens_id


def extract_magma_z(magma_z_score_file, gene_id_to_ens_id,trait_name):
	f = open(magma_z_score_file)
	ens_id_to_magma_z = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			matched_indices = np.where(np.asarray(data)==trait_name)[0]
			if len(matched_indices) != 1:
				print('assumptino error')
				pdb.set_trace()
			trait_index = matched_indices[0]
			continue
		gene_name = data[0]
		magma_z = data[trait_index]

		if gene_name not in gene_id_to_ens_id:
			continue
		ens_id = gene_id_to_ens_id[gene_name]


		if ens_id in ens_id_to_magma_z:
			print('asssumptioenroroo')
			pdb.set_trace()

		if magma_z =='NA':
			ens_id_to_magma_z[ens_id] = np.nan
		else:
			ens_id_to_magma_z[ens_id] = float(magma_z)

	f.close()

	return ens_id_to_magma_z


def extract_sgd_links_data(sgdlinks_gene_summary_file, filter_means_thresh=None):
	genes = []
	zs = []
	means = []
	used_genes = {}

	f = open(sgdlinks_gene_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ens_id = data[0].split('.')[0]
		score = float(data[1])
		z_score = float(data[5])

		if ens_id in used_genes:
			print('assumption eroror')
			pdb.set_trace()
		used_genes[ens_id] = 1

		if filter_means_thresh is not None and score > filter_means_thresh:
			continue

		genes.append(ens_id)
		zs.append(z_score)
		means.append(score)

	f.close()

	return np.asarray(genes), np.asarray(zs), np.asarray(means)


def run_pops_enrichment_meta_analysis(gene_set_enrichment_results_dir, method_identifier, gwas_trait_names, meta_analysis_output_stem, data_type):

	n_genes_arr = [5,10,20,50,100,200,500]

	sgds_dicti = {}
	magmas_dicti = {}

	for n_genes in n_genes_arr:
		sgds_dicti[n_genes] = []
		magmas_dicti[n_genes] = []

	for trait_name in gwas_trait_names:
		genes = []
		sgds = []
		magmas = []
		pops = []

		enrichment_summary_file = gene_set_enrichment_results_dir + trait_name + '_' + method_identifier + '_' + data_type + '_pops_enrichments_summary.txt'
		f = open(enrichment_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			genes.append(data[0])
			sgds.append(float(data[1]))
			magmas.append(float(data[2]))
			pops.append(float(data[3]))
		f.close()

		# Organize data
		genes = np.asarray(genes)
		sgds = np.asarray(sgds)
		magmas = np.asarray(magmas)
		pops = np.asarray(pops)

		for n_genes in n_genes_arr:

			pops_sgds_arr = pops[np.argsort(-sgds)][:n_genes]
			sgds_dicti[n_genes].append(pops_sgds_arr)
			pops_magmas_arr = pops[np.argsort(-magmas)][:n_genes]
			magmas_dicti[n_genes].append(pops_magmas_arr)


	pops_mean_output_file = meta_analysis_output_stem + '_' + data_type + '_pops_enrichments_avg_pops_scores_per_threshold.txt'
	t = open(pops_mean_output_file,'w')
	t.write('method\tn_genes\tPOPS_mean\tSE\tPOPS_lb\tPOPS_ub\n')

	for n_genes in n_genes_arr:

		pops_arr = np.hstack(sgds_dicti[n_genes])
		mean = np.mean(pops_arr)
		se = np.std(pops_arr)/np.sqrt(len(pops_arr))
		lb = mean - (1.96*se)
		ub = mean + (1.96*se)

		t.write('sgdLinks\t' + str(n_genes) + '\t' + str(mean) + '\t' + str(se) + '\t' + str(lb) + '\t' + str(ub) + '\n')

		pops_arr = np.hstack(magmas_dicti[n_genes])
		mean = np.mean(pops_arr)
		se = np.std(pops_arr)/np.sqrt(len(pops_arr))
		lb = mean - (1.96*se)
		ub = mean + (1.96*se)

		t.write('Magma\t' + str(n_genes) + '\t' + str(mean) + '\t' + str(se) + '\t' + str(lb) + '\t' + str(ub) + '\n')

	t.close()

	print(pops_mean_output_file)

	return



def run_pops_enrichment_meta_analysis_z_thresh(gene_set_enrichment_results_dir, method_identifier, gwas_trait_names, meta_analysis_output_stem, data_type):

	genes = []
	sgds = []
	magmas = []
	pops = []

	for trait_name in gwas_trait_names:


		enrichment_summary_file = gene_set_enrichment_results_dir + trait_name + '_' + method_identifier + '_' + data_type + '_pops_enrichments_summary.txt'
		f = open(enrichment_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			genes.append(data[0])
			sgds.append(float(data[1]))
			magmas.append(float(data[2]))
			pops.append(float(data[3]))
		f.close()

	# Organize data
	genes = np.asarray(genes)
	sgds = np.asarray(sgds)
	magmas = np.asarray(magmas)
	pops = np.asarray(pops)


	pops_mean_output_file = meta_analysis_output_stem + '_avg_pops_scores_per_z_score_threshold.txt'
	t = open(pops_mean_output_file,'w')
	t.write('method\tz_threshold\tn_genes\tPOPS_mean\tSE\tPOPS_lb\tPOPS_ub\n')

	z_thresholds = [6.0, 5.0, 4.0, 3.0, 2.5, 2.0, 1.5, 1.0]

	for z_thresh in z_thresholds:

		indices = sgds > z_thresh
		n_genes = np.sum(indices)
		pops_arr = pops[indices]
		mean = np.mean(pops_arr)
		se = np.std(pops_arr)/np.sqrt(len(pops_arr))
		lb = mean - (1.96*se)
		ub = mean + (1.96*se)
		t.write('sgdLinks\t' + str(z_thresh) + '\t' + str(n_genes) + '\t' + str(mean) + '\t' + str(se) + '\t' + str(lb) + '\t' + str(ub) + '\n')

		pops_arr = pops[np.argsort(-magmas)][:n_genes]
		mean = np.mean(pops_arr)
		se = np.std(pops_arr)/np.sqrt(len(pops_arr))
		lb = mean - (1.96*se)
		ub = mean + (1.96*se)

		t.write('Magma\t' + str(z_thresh) + '\t' + str(n_genes) + '\t' + str(mean) + '\t' + str(se) + '\t' + str(lb) + '\t' + str(ub) + '\n')

	t.close()

	print(pops_mean_output_file)

	return



def run_pops_enrichment_analysis(sgdlinks_genes, sgdlinks_scores, ens_id_to_magma_z, ens_id_to_pops_score, sgdlink_data_name, output_stem):
	# Keep track of data
	genes = []
	sgds = []
	magmas = []
	pops = []

	# First create enrichment_summary file
	enrichment_summary_file = output_stem + '_summary.txt'
	t = open(enrichment_summary_file,'w')
	t.write('ens_id\t' + sgdlink_data_name + '\t' + 'MAGMA_z\tPOPS_score\n')

	for ii, ens_id in enumerate(sgdlinks_genes):
		sgdlinks_score = sgdlinks_scores[ii]
		if ens_id not in ens_id_to_pops_score or ens_id not in ens_id_to_magma_z:
			continue
		pops_score = ens_id_to_pops_score[ens_id]
		magma_z = ens_id_to_magma_z[ens_id]
		if np.isnan(pops_score) or np.isnan(magma_z):
			continue
		t.write(ens_id + '\t' + str(sgdlinks_score) + '\t' + str(magma_z) + '\t' + str(pops_score) + '\n')
		genes.append(ens_id)
		sgds.append(sgdlinks_score)
		magmas.append(magma_z)
		pops.append(pops_score)
	t.close()

	# Organize data
	genes = np.asarray(genes)
	sgds = np.asarray(sgds)
	magmas = np.asarray(magmas)
	pops = np.asarray(pops)


	pops_mean_output_file = output_stem + '_avg_pops_scores_per_threshold.txt'
	t = open(pops_mean_output_file,'w')
	t.write('method\tn_genes\tPOPS_mean\tSE\tPOPS_lb\tPOPS_ub\n')

	for n_genes in [5,10,20,50,100,200,500]:

		pops_arr = pops[np.argsort(-sgds)][:n_genes]
		mean = np.mean(pops_arr)
		se = np.std(pops_arr)/np.sqrt(len(pops_arr))
		lb = mean - (1.96*se)
		ub = mean + (1.96*se)

		t.write('sgdLinks\t' + str(n_genes) + '\t' + str(mean) + '\t' + str(se) + '\t' + str(lb) + '\t' + str(ub) + '\n')

		pops_arr = pops[np.argsort(-magmas)][:n_genes]
		mean = np.mean(pops_arr)
		se = np.std(pops_arr)/np.sqrt(len(pops_arr))
		lb = mean - (1.96*se)
		ub = mean + (1.96*se)

		t.write('Magma\t' + str(n_genes) + '\t' + str(mean) + '\t' + str(se) + '\t' + str(lb) + '\t' + str(ub) + '\n')

	t.close()

	return

def extract_ldl_silverstandard_genes(ldl_silverstandard_gene_set, gene_id_to_ens_id):
	dicti = {}
	f = open(ldl_silverstandard_gene_set)
	tmp_mapping = {}
	tmp_mapping['known'] = 1
	tmp_mapping['bystander'] = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if len(data) != 12:
			print('assumtpioneroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		if gene_id not in gene_id_to_ens_id:
			continue
		ens_id = gene_id_to_ens_id[gene_id]
		gene_type = data[11]
		if ens_id in dicti:
			print('assumption eroorro')
			pdb.set_trace()
		dicti[ens_id] = tmp_mapping[gene_type]

	f.close()




	return dicti

def run_ldl_silver_standard_geneset_enrichment_analysis(sgdlinks_genes, sgdlinks_scores, ens_id_to_magma_z, ens_id_to_pops_score, ens_id_to_ldl_gs, sgdlink_data_name, output_stem):
	# Keep track of data
	genes = []
	sgds = []
	magmas = []
	pops = []
	ldl_gs = []

	# First create enrichment_summary file
	enrichment_summary_file = output_stem + '_summary.txt'
	t = open(enrichment_summary_file,'w')
	t.write('ens_id\t' + sgdlink_data_name + '\t' + 'MAGMA_z\tPOPS_score\tldl_silverstandard_gs\n')

	for ii, ens_id in enumerate(sgdlinks_genes):
		sgdlinks_score = sgdlinks_scores[ii]
		if ens_id not in ens_id_to_pops_score or ens_id not in ens_id_to_magma_z or ens_id not in ens_id_to_ldl_gs:
			continue
		pops_score = ens_id_to_pops_score[ens_id]
		magma_z = ens_id_to_magma_z[ens_id]
		ldl_g = ens_id_to_ldl_gs[ens_id]
		if np.isnan(pops_score) or np.isnan(magma_z):
			continue
		t.write(ens_id + '\t' + str(sgdlinks_score) + '\t' + str(magma_z) + '\t' + str(pops_score) + '\t' + str(ldl_g) + '\n')
		genes.append(ens_id)
		sgds.append(sgdlinks_score)
		magmas.append(magma_z)
		pops.append(pops_score)
		ldl_gs.append(ldl_g)
	t.close()

	# Make FDR-power curve
	fdr_power_ouptut_file = output_stem + '_fdr_power_curve.txt'
	t = open(fdr_power_ouptut_file,'w')
	t.write('method\tthreshold\tfdr\tpower\n')

	# Organize data
	genes = np.asarray(genes)
	sgds = np.asarray(sgds)
	magmas = np.asarray(magmas)
	pops = np.asarray(pops)
	ldl_gs = np.asarray(ldl_gs)

	total_pos = np.sum(ldl_gs == 1)

	for sgd_thresh in np.sort(sgds):
		indices = sgds >= sgd_thresh
		vals = ldl_gs[indices]
		fdr = np.sum(vals==0)/len(vals)
		power = np.sum(vals==1)/total_pos
		t.write('sgdLinks\t' + str(sgd_thresh) + '\t' + str(fdr) + '\t' + str(power) + '\n')

	for magma_thresh in np.sort(magmas):
		indices = magmas >= magma_thresh
		vals = ldl_gs[indices]
		fdr = np.sum(vals==0)/len(vals)
		power = np.sum(vals==1)/total_pos
		t.write('Magma\t' + str(magma_thresh) + '\t' + str(fdr) + '\t' + str(power) + '\n')	

	t.close()

	print(fdr_power_ouptut_file)
	'''
	for n_genes in [1, 2, 4, 6, 10, 15, 20, 40, 100, 200, len(genes)]:
		print(str(n_genes) + '\tsgd\t' + str(np.mean(ldl_gs[np.argsort(-sgds)][:n_genes])) + '\t' + str(np.sum(ldl_gs[np.argsort(-sgds)][:n_genes])))
		print(str(n_genes) + '\tmagma\t' + str(np.mean(ldl_gs[np.argsort(-magmas)][:n_genes])) + '\t' + str(np.sum(ldl_gs[np.argsort(-magmas)][:n_genes])))
	'''
	return






def extract_gwas_trait_names(gwas_traits_file):
	f = open(gwas_traits_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1 
			continue
		arr.append(data[0])
	f.close()

	return np.asarray(arr)




######################
# Command line args
######################
gwas_traits_file = sys.argv[1]
learned_snp_gene_links_dir = sys.argv[2]
pops_results_summary_file = sys.argv[3]
magma_z_score_file = sys.argv[4]
ldl_silverstandard_gene_set = sys.argv[5]
gene_set_enrichment_results_dir = sys.argv[6]
method_identifier = sys.argv[7]




# Extract gwas trait names
gwas_trait_names = extract_gwas_trait_names(gwas_traits_file)






#########################
# POPS gene set enrichment analysis
#########################
# Loop through gwas traits
for trait_name in gwas_trait_names:

	# Per trait output stem
	per_trait_output_stem = gene_set_enrichment_results_dir + trait_name + '_' + method_identifier

	# Create mapping from trait name to pops trait name
	pops_trait_mapping = create_dictionary_mapping_trait_names_to_pops_trait_names()
	pops_trait_name = pops_trait_mapping[trait_name]


	# Extract pops scores
	ens_id_to_pops_score, gene_id_to_ens_id = extract_pops_scores(pops_results_summary_file, pops_trait_name)

	# Extract magma z-scores
	ens_id_to_magma_z = extract_magma_z(magma_z_score_file, gene_id_to_ens_id,trait_name)

	# Extract sgdlinks genes and scores
	sgdlinks_gene_summary_file = learned_snp_gene_links_dir + 'snp_gene_links_' + trait_name + '_lmm_snp_gene_link_' + method_identifier + '_gene_score_averaged.txt'
	sgdlinks_genes, sgdlinks_zs, sgdlinks_avg_vars = extract_sgd_links_data(sgdlinks_gene_summary_file, filter_means_thresh=1.0)


	# Run enrichment analysis with respect to sgdlinks zs
	output_stem = per_trait_output_stem + '_zscore_pops_enrichments'
	run_pops_enrichment_analysis(sgdlinks_genes, sgdlinks_zs, ens_id_to_magma_z, ens_id_to_pops_score, 'sgdLinks_z', output_stem)

# Run POPS enrichment meta-analysis
meta_analysis_output_stem = gene_set_enrichment_results_dir + 'meta_analysis_' + method_identifier
run_pops_enrichment_meta_analysis(gene_set_enrichment_results_dir, method_identifier, gwas_trait_names, meta_analysis_output_stem, 'zscore')


# Run POPS enrichment meta-analysis
meta_analysis_output_stem = gene_set_enrichment_results_dir + 'meta_analysis_' + method_identifier
run_pops_enrichment_meta_analysis_z_thresh(gene_set_enrichment_results_dir, method_identifier, gwas_trait_names, meta_analysis_output_stem, 'zscore')





#########################
# Silver-standard LDL cholesterol gene set enrichment analysis
#########################
trait_name = 'UKB_460K.biochemistry_LDLdirect'
per_trait_output_stem = gene_set_enrichment_results_dir + trait_name + '_' + method_identifier

# Create mapping from trait name to pops trait name
pops_trait_mapping = create_dictionary_mapping_trait_names_to_pops_trait_names()
pops_trait_name = pops_trait_mapping[trait_name]


# Extract pops scores
ens_id_to_pops_score, gene_id_to_ens_id = extract_pops_scores(pops_results_summary_file, pops_trait_name)

# Extract magma z-scores
ens_id_to_magma_z = extract_magma_z(magma_z_score_file, gene_id_to_ens_id,trait_name)

# Extract sgdlinks genes and scores
sgdlinks_gene_summary_file = learned_snp_gene_links_dir + 'snp_gene_links_' + trait_name + '_lmm_snp_gene_link_' + method_identifier + '_gene_score_averaged.txt'
sgdlinks_genes, sgdlinks_zs, sgdlinks_avg_vars = extract_sgd_links_data(sgdlinks_gene_summary_file,filter_means_thresh=1.0)

# Extract LDL silverstandard geneset
ens_id_to_ldl_silverstandard = extract_ldl_silverstandard_genes(ldl_silverstandard_gene_set, gene_id_to_ens_id)


# Run enrichment analysis with respect to sgdlinks zs
output_stem = per_trait_output_stem + '_zscore_ldl_silver_standard_gs_enrichments'
run_ldl_silver_standard_geneset_enrichment_analysis(sgdlinks_genes, sgdlinks_zs, ens_id_to_magma_z, ens_id_to_pops_score, ens_id_to_ldl_silverstandard, 'sgdLinks_z', output_stem)

