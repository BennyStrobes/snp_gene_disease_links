import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np
import os
import pdb
import gzip
from pandas_plink import read_plink1_bin




def extract_dictionary_list_of_hapmap3_rsids(hapmap3_rsid_file):
	f = open(hapmap3_rsid_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] in dicti:
			print('assumptioneroeor')
			pdb.set_trace()
		dicti[data[0]] = 1
	f.close()
	return dicti


def extract_dictionary_list_of_reference_rsids(filer):
	f = gzip.open(filer)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[2]
		dicti[rsid] = 1
	f.close()
	return dicti


def save_snp_rsids_to_text_file(snp_id_output_file, window_rsids):
	t = open(snp_id_output_file,'w')
	t.write('rsid\n')
	for window_rsid in window_rsids:
		t.write(window_rsid + '\n')
	t.close()
	return

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


######################
# Command line args
######################
chrom_num = sys.argv[1]
genome_wide_ld_windows_file = sys.argv[2]
hapmap3_rsid_file = sys.argv[3]
baselineLD_anno_dir = sys.argv[4]
kg_plink_dir = sys.argv[5]
snp_set = sys.argv[6]
output_dir = sys.argv[7]


# Contains hapmap3 rsids from all chromosomes (this should be cool)
regression_rsid_dictionary = extract_dictionary_list_of_hapmap3_rsids(hapmap3_rsid_file)

# Contains all reference rsids present in 1K genomes
reference_rsid_dictionary = extract_dictionary_list_of_reference_rsids(baselineLD_anno_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz')


# Load in Reference Genotype data
geno_stem = kg_plink_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.'

G_obj = read_plink1_bin(geno_stem + 'bed', geno_stem + 'bim', geno_stem + 'fam', verbose=False)
G = G_obj.values # Numpy 2d array of dimension num samples X num snps
ref_chrom = np.asarray(G_obj.chrom)
ref_pos = np.asarray(G_obj.pos)
# For our purposes, a0 is the effect allele
# For case of plink package, a0 is the first column in the plink bim file
ref_a0 = np.asarray(G_obj.a0)
ref_a1 = np.asarray(G_obj.a1)
snp_cm = np.asarray(G_obj.cm)
snp_rs_ids = np.asarray(G_obj.snp)
n_ref_samples = G.shape[0]

# Filter reference data to snps that are reference snps
valid_snps = []
for rsid in snp_rs_ids:
	if rsid in reference_rsid_dictionary:
		valid_snps.append(True)
	else:
		valid_snps.append(False)
valid_snps = np.asarray(valid_snps)
if len(valid_snps) != np.sum(valid_snps):
	print('assumption eroror')
	pdb.set_trace()


# Loop through windows
# For each window extract LD
f = open(genome_wide_ld_windows_file)
# Also open output file
output_file = output_dir + 'window_LD_summary_' + snp_set + '_chr' + str(chrom_num) + '.txt'
t = open(output_file,'w')
t.write('window_id\tchrom_num\tstart_position\tend_position\tLD_mat\tsnp_file\tLD_EIVD_Q_mat_file\tLD_EIVD_W_mat_file\n')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Skip windows not on this chromosome
	line_chrom_num = data[0]
	if line_chrom_num != str(chrom_num):
		continue

	# Extract relevent fields
	window_start_pos = int(data[1])
	window_end_pos = int(data[2])

	window_name = line_chrom_num + ':' + str(window_start_pos) + ':' + str(window_end_pos)
	print(window_name)

	# Extract window indices
	window_indices = (ref_pos >= window_start_pos) & (ref_pos < window_end_pos)

	# Filter to window
	window_rsids = snp_rs_ids[window_indices]
	window_snp_pos = ref_pos[window_indices]
	window_G = G[:, window_indices]

	# Compute window LD
	window_LD = np.corrcoef(np.transpose(window_G))
	print(window_LD.shape)

	# Extract indice of snps that are regression snps
	if snp_set == 'hampmap3_snps':
		regression_snp_boolean = []
		for ii, window_rsid in enumerate(window_rsids):
			tmp_pos = window_snp_pos[ii]
			if window_rsid in regression_rsid_dictionary:
				regression_snp_boolean.append(True)
			else:
				regression_snp_boolean.append(False)
		regression_snp_boolean = np.asarray(regression_snp_boolean)
	else:
		print('assumption erororo')
		print('not yet implemented')

	# Skip windows with very few regression snps
	if np.sum(regression_snp_boolean) < 5:
		print('skip window')
		continue

	# Filter snps to only regression snps
	window_final_rsids = window_rsids[regression_snp_boolean]
	filtered_LD = window_LD[regression_snp_boolean, :][:, regression_snp_boolean]


	# Run eigenvalue decomposition on LD matrix
	Q_mat, w_premult = eigenvalue_decomp_ld(filtered_LD)
	print(Q_mat.shape)

	# Save LD file to npy file
	ld_output_file = output_dir + snp_set + '_' + window_name + '_LD.npy'
	np.save(ld_output_file, filtered_LD)

	# Save snp rsids to text file
	snp_id_output_file = output_dir + snp_set + '_' + window_name + '_rsids.txt'
	save_snp_rsids_to_text_file(snp_id_output_file, window_final_rsids)

	# Save EIVD Q mat file to npy file
	Q_output_file = output_dir + snp_set + '_' + window_name + '_LD_EIVD_Q_mat.npy'
	np.save(Q_output_file, Q_mat)
	
	# Save EIVD w mat file to npy file
	W_output_file = output_dir + snp_set + '_' + window_name + '_LD_EIVD_W_mat.npy'
	np.save(W_output_file, w_premult)

	t.write(window_name + '\t' + str(chrom_num) + '\t' + str(window_start_pos) + '\t' + str(window_end_pos) + '\t' + ld_output_file + '\t' + snp_id_output_file + '\t' + Q_output_file + '\t' + W_output_file + '\n')

f.close()
t.close()






