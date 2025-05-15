import numpy as np
import os
import sys
import pdb


def extract_rsids_and_alleles(snp_file):
	rsids = []
	a0s = []
	a1s = []
	f = open(snp_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsids.append(data[0])
		a0s.append(data[1])
		a1s.append(data[2])

	f.close()


	return np.asarray(rsids), np.asarray(a0s), np.asarray(a1s)




def load_in_rsid_to_alleles(input_window_summary_file):
	rsid_to_alleles = {}
	head_count = 0
	f = open(input_window_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_file = data[5]
		rsids, a0s, a1s = extract_rsids_and_alleles(snp_file)

		for ii, rsid in enumerate(rsids):
			if rsid in rsid_to_alleles:
				print('assumpriont oeororor')
				pdb.set_trace()
			rsid_to_alleles[rsid] = [a0s[ii], a1s[ii]]
	f.close()
	return rsid_to_alleles




def create_dictionary_mapping_from_rsid_to_z(trait_sumstat_file, rsid_to_alleles):
	f = open(trait_sumstat_file)
	rsid_to_z = {}
	skipped = 0
	counter = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_rsid = data[0]
		line_a0 = data[1]
		line_a1 = data[2]
		line_z = float(data[4])
		gwas_ss = float(data[3])

		if line_rsid not in rsid_to_alleles:
			continue

		ref_a0 = rsid_to_alleles[line_rsid][0]
		ref_a1 = rsid_to_alleles[line_rsid][1]

		if ref_a0 == line_a0 and ref_a1 == line_a1:
			line_z = line_z*1.0
		elif ref_a0 == line_a1 and ref_a1 == line_a0:
			line_z = line_z*-1.0
		else:
			skipped = skipped +1
			continue
		if line_rsid in rsid_to_z:
			print('assumption oeroror')
			pdb.set_trace()
		rsid_to_z[line_rsid] = line_z

	f.close()


	return rsid_to_z, gwas_ss


def get_z_scores_for_this_window(rsid_file, rsid_to_z):
	f = open(rsid_file)
	window_zs = []
	window_snp_filter = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_rsid = data[0]
		if line_rsid not in rsid_to_z:
			window_snp_filter.append(False)
		else:
			window_snp_filter.append(True)
			window_zs.append(rsid_to_z[line_rsid])
	f.close()

	return np.asarray(window_zs), np.asarray(window_snp_filter)



def compute_lambda_thresh(lambdas, rho_thresh):
	totaler = np.sum(lambdas)
	cur_total = 0
	lambda_thresh = -1
	for lambda_val in -np.sort(-lambdas):
		cur_total = cur_total + lambda_val
		if cur_total/totaler > rho_thresh:
			if lambda_thresh == -1:
				lambda_thresh = lambda_val


	if lambda_thresh == -1:
		print('assumption eroror')
		pdb.set_trace()

	return lambda_thresh


def eigenvalue_decomp_ld(ld_mat):
	# EIG value decomp
	lambdas_full, U_full = np.linalg.eig(ld_mat)
	non_negative_components = lambdas_full > 0.0
	lambdas = lambdas_full[non_negative_components]
	U = U_full[:, non_negative_components]
	real_components = np.iscomplex(lambdas) == False
	lambdas = lambdas[real_components]
	U = U[:, real_components]
	if np.sum(np.iscomplex(lambdas)) > 0:
		print('assumption eroror')
		pdb.set_trace()
	lambdas = lambdas.astype(float)
	U = U.astype(float)

	rho_thresh = 0.99
	lambda_thresh = compute_lambda_thresh(lambdas, rho_thresh)
	thresh_components = lambdas >= lambda_thresh
	lambdas = lambdas[thresh_components]
	U = U[:, thresh_components]


	# Note that reconstruction of ld_mat is achieved with np.dot(np.dot(U, np.diag(lambdas)), np.transpose(U))

	# Compute some relevent quantities
	Q_mat = np.dot(np.diag(lambdas**(.5)), np.transpose(U))
	w_premult = np.dot(np.diag(lambdas**(-.5)), np.transpose(U))    

	return Q_mat, w_premult


def create_filtered_rsid_file(rsid_file, new_rsid_file, window_snp_filter):
	f = open(rsid_file)
	t = open(new_rsid_file,'w')
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		if window_snp_filter[counter]:
			t.write(line + '\n')
		counter = counter + 1

	if counter != len(window_snp_filter):
		print('assumption eroror')
		pdb.set_trace()
	f.close()
	t.close()
	return



#######################
# Command line args
#######################
trait_sumstat_file = sys.argv[1]
input_window_summary_file = sys.argv[2]
gene_summary_file = sys.argv[3]
disease_specific_tmp_data_dir = sys.argv[4]


######
# Load in rsid to alleles
rsid_to_alleles = load_in_rsid_to_alleles(input_window_summary_file)


######
# Load in sumstats (get alleles in line with LD!)
# Create dictionary mapping from rsid to z-scores
# Also extract gwas sample size
rsid_to_z, gwas_sample_size = create_dictionary_mapping_from_rsid_to_z(trait_sumstat_file, rsid_to_alleles)


t = open(disease_specific_tmp_data_dir + '_window_summary_file.txt','w')

# Loop through windows
# For each window, reprocesss data

f = open(input_window_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + 'gwas_sample_size\t')
		t.write(data[5] + '\t' + 'EIVD_gwas_z_file' + '\t' + data[6] + '\t' + data[8] + '\t' + data[9] + '\n')
		continue

	# Extract relevent fields
	window_name = data[0]
	print(window_name)
	LD_file = data[4]
	rsid_file = data[5]
	EIVD_Q_file = data[6]
	EIVD_W_file = data[7]
	linked_genes_file = data[8]
	snp_gene_anno_file = data[9]


	# Get z-scores for this window
	window_z, window_snp_filter = get_z_scores_for_this_window(rsid_file, rsid_to_z)
	if len(window_z) < 5:
		continue

	# Filter LD matrix to have correct fields
	orig_LD = np.load(LD_file)
	new_LD = orig_LD[:, window_snp_filter][window_snp_filter, :]

	# EIVD on new LD matrix
	Q_mat, w_premult = eigenvalue_decomp_ld(new_LD)

	# Get gwas z-scores in transformed space
	EIVD_zeds = np.dot(w_premult, window_z)
	new_EIVD_zed_file = disease_specific_tmp_data_dir + '_' + window_name + '_EIVD_zeds.npy'
	np.save(new_EIVD_zed_file, EIVD_zeds)

	# Create new, filtered rsid file
	new_rsid_file = disease_specific_tmp_data_dir + '_' + window_name + '_rsids.txt'
	create_filtered_rsid_file(rsid_file, new_rsid_file, window_snp_filter)

	# Save new Q mat file
	new_EIVD_Q_file = disease_specific_tmp_data_dir + '_' + window_name + '_LD_EIVD_Q_mat.npy'
	np.save(new_EIVD_Q_file, Q_mat)

	# Load in linked genes data
	linked_genes = np.load(linked_genes_file)
	# Filter linked genes
	new_linked_genes = linked_genes[window_snp_filter, :]
	# Save new linked genes
	new_linked_genes_file = disease_specific_tmp_data_dir + '_' + window_name + '_linked_integer_genes.npy'
	np.save(new_linked_genes_file, new_linked_genes)

	# Load in snp_gene_anno
	snp_gene_anno = np.load(snp_gene_anno_file)
	new_snp_gene_anno = snp_gene_anno[window_snp_filter, :, :]
	# Save new snp gene anno
	new_snp_gene_anno_file = disease_specific_tmp_data_dir + '_' + window_name + '_snp_gene_anno.npy'
	np.save(new_snp_gene_anno_file, new_snp_gene_anno)


	# Write to output file
	t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + str(gwas_sample_size) + '\t')
	t.write(new_rsid_file + '\t' + new_EIVD_zed_file + '\t' + new_EIVD_Q_file + '\t' + new_linked_genes_file + '\t' + new_snp_gene_anno_file + '\n')
f.close()
t.close()





