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
	"""Map each regression snp to its reference alleles, and record each window's snp order.

	Column 5 (regression_snp_file) is deliberate. Gwas z-scores attach to the regression
	(output) snps, which index the ROWS of Q_mat, so those are the snps needing allele
	alignment. Column 4 is the much larger reference (input) snp set indexing Q_mat's
	COLUMNS; it carries the gene annotations and never gets z-scores. Do not "correct"
	this to data[4] -- preorganize_snp_gene_annotations.py needs the opposite column,
	for exactly the same reason.

	window_to_rsids is returned so the per-window z-score lookup does not have to re-read
	every window's rsid file a second time.
	"""
	rsid_to_alleles = {}
	window_to_rsids = {}
	head_count = 0
	f = open(input_window_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# 8 LD-summary columns plus linked_genes_file and snp_gene_anno_file
		if len(data) != 10:
			print('assumption eroror: expected 10 columns in ' + input_window_summary_file + ', got ' + str(len(data)))
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		regression_snp_file = data[5]
		rsids, a0s, a1s = extract_rsids_and_alleles(regression_snp_file)

		# Canonical z-score ordering for this window (rows of Q_mat)
		window_to_rsids[window_name] = rsids

		for ii in range(len(rsids)):
			rsid = str(rsids[ii])
			if rsid in rsid_to_alleles:
				print('assumption eroror: regression snp ' + rsid + ' appears in more than one window')
				pdb.set_trace()
			# plain str/tuple rather than numpy scalars in a list: measured 0.36 GB -> 0.13 GB
			# across a ~1.2M snp regression set
			rsid_to_alleles[rsid] = (str(a0s[ii]), str(a1s[ii]))
	f.close()
	return rsid_to_alleles, window_to_rsids




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


def get_z_scores_for_this_window(window_rsids, rsid_to_z):
	"""Z-scores for one window, in the canonical regression-snp order.

	Returns a FULL LENGTH vector (one entry per regression snp, ie. per column of the
	window's whitener) with np.nan wherever this trait has no usable z-score, plus a
	boolean mask of which entries were actually observed. Keeping the vector full length
	is what lets every trait keep sharing the same w_premult / Q_mat.
	"""
	window_zs = []
	observed_boolean = []
	for rsid in window_rsids:
		rsid = str(rsid)
		if rsid not in rsid_to_z:
			window_zs.append(np.nan)
			observed_boolean.append(False)
		else:
			window_zs.append(rsid_to_z[rsid])
			observed_boolean.append(True)

	return np.asarray(window_zs), np.asarray(observed_boolean)


def reconstruct_regression_ld_from_whitener(w_premult):
	"""Rebuild the regression-snp LD matrix R_OO from the saved whitener.

	generate_ld_matrices_for_ld_windows_on_single_chromosome.py deliberately stops saving
	R_OO, but it does save w_premult = Lambda^(-.5) U'. The columns of its transpose are
	U_i scaled by lambda_i^(-.5), so both U and Lambda drop straight back out and
	R_OO = U Lambda U'.

	This is the rank-k reconstruction left after that script's rho=0.99 truncation, ie. an
	already-smoothed R_OO. Harmless (helpful, even) for imputation, which regularizes
	anyway -- but it is why the diagonal is not exactly 1.
	"""
	w_transpose = np.transpose(w_premult)
	inv_sqrt_lambdas = np.sqrt(np.sum(np.square(w_transpose), axis=0))
	if np.sum(inv_sqrt_lambdas <= 0.0) > 0:
		print('assumption eroror: non-positive eigenvalue recovered from whitener')
		pdb.set_trace()
	lambdas = 1.0/np.square(inv_sqrt_lambdas)
	U = w_transpose/inv_sqrt_lambdas

	return np.dot(U*lambdas, np.transpose(U))


def impute_missing_z(window_z, observed_boolean, regression_ld, ridge_lambda):
	"""Ridge-regularized conditional-mean imputation of missing z-scores.

		z_miss = R_mo (R_oo + lambda I)^-1 z_obs

	The ridge is not optional. R_oo is estimated from a ~500 sample reference panel while
	windows carry far more regression snps than that, so R_oo is rank deficient and has no
	plain inverse.

	Also returns the per-imputed-snp imputation r2 (the diagonal of R_mo (R_oo+lambda I)^-1 R_om).
	Near 1 means the snp is tightly tagged by its observed neighbours and the imputed value is
	trustworthy; near 0 means it is effectively independent of everything observed, and the
	imputed value is close to made up.
	"""
	missing_boolean = observed_boolean == False

	R_oo = regression_ld[observed_boolean, :][:, observed_boolean]
	R_mo = regression_ld[missing_boolean, :][:, observed_boolean]
	R_oo_ridge = R_oo + ridge_lambda*np.eye(R_oo.shape[0])

	# One factorization each: R_oo^-1 applied to the observed z, and to R_om
	solved_z = np.linalg.solve(R_oo_ridge, window_z[observed_boolean])
	solved_R = np.linalg.solve(R_oo_ridge, np.transpose(R_mo))

	imputation_r2 = np.sum(np.transpose(solved_R)*R_mo, axis=1)

	completed_z = np.copy(window_z)
	completed_z[missing_boolean] = np.dot(R_mo, solved_z)

	return completed_z, imputation_r2



#######################
# Command line args
#######################
trait_sumstat_file = sys.argv[1]
input_window_summary_file = sys.argv[2]
gene_summary_file = sys.argv[3]
disease_specific_tmp_data_dir = sys.argv[4]


# Ridge added to the diagonal of R_oo before inverting (the ImpG-Summary convention).
ridge_lambda = 0.1

# Windows retaining too little of the regression set are skipped rather than imputed:
# predicting most of a window from a small remainder is not meaningful. With a correctly
# built all-trait intersection these should never fire.
min_observed_fraction = 0.5
min_observed_snps = 5


######
# Load in rsid to alleles (and the canonical per-window regression snp ordering)
rsid_to_alleles, window_to_rsids = load_in_rsid_to_alleles(input_window_summary_file)


######
# Load in sumstats (get alleles in line with LD!)
# Create dictionary mapping from rsid to z-scores
# Also extract gwas sample size
rsid_to_z, gwas_sample_size = create_dictionary_mapping_from_rsid_to_z(trait_sumstat_file, rsid_to_alleles)


t = open(disease_specific_tmp_data_dir + '_window_summary_file.txt','w')

# Loop through windows.
# Everything except the z-scores is trait invariant, so this step writes exactly ONE array
# per window (the whitened z vector) and REFERENCES the shared artifacts for everything else.

f = open(input_window_summary_file)
head_count = 0

n_windows_used = 0
n_windows_skipped = 0
n_snps_total = 0
n_snps_imputed = 0
imputation_r2_arr = []

for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + 'gwas_sample_size\t')
		t.write(data[4] + '\t' + 'EIVD_gwas_z_file' + '\t' + data[6] + '\t' + data[8] + '\t' + data[9] + '\n')
		continue

	# Shared artifacts are REFERENCED, never copied:
	#   data[4] reference rsids -> index the COLUMNS of Q_mat (ie. gamma, and n_snps downstream)
	#   data[6] Q_mat  data[7] w_premult  data[8] linked genes  data[9] snp-gene anno
	window_name = data[0]
	reference_rsid_file = data[4]
	EIVD_Q_file = data[6]
	EIVD_W_file = data[7]
	linked_genes_file = data[8]
	snp_gene_anno_file = data[9]

	window_rsids = window_to_rsids[window_name]

	# Get z-scores for this window (np.nan where this trait has no usable z)
	window_z, observed_boolean = get_z_scores_for_this_window(window_rsids, rsid_to_z)

	n_observed = np.sum(observed_boolean)
	n_missing = len(observed_boolean) - n_observed

	if n_observed < min_observed_snps or n_observed < min_observed_fraction*len(observed_boolean):
		print('skipping window ' + window_name + ': only ' + str(n_observed) + ' of ' + str(len(observed_boolean)) + ' regression snps observed')
		n_windows_skipped = n_windows_skipped + 1
		continue

	w_premult = np.load(EIVD_W_file)

	if n_missing > 0:
		# Ridge conditional-mean imputation, off R_OO rebuilt from the saved whitener
		regression_ld = reconstruct_regression_ld_from_whitener(w_premult)
		window_z, imputation_r2 = impute_missing_z(window_z, observed_boolean, regression_ld, ridge_lambda)
		imputation_r2_arr.append(imputation_r2)
		n_snps_imputed = n_snps_imputed + n_missing

	if np.sum(np.isnan(window_z)) > 0:
		print('assumption eroror: window ' + window_name + ' still carries missing z-scores after imputation')
		pdb.set_trace()

	n_snps_total = n_snps_total + len(window_z)
	n_windows_used = n_windows_used + 1

	# Whiten. w_premult is trait invariant, so this is the ONLY per-trait array produced.
	EIVD_zeds = np.dot(w_premult, window_z)
	new_EIVD_zed_file = disease_specific_tmp_data_dir + '_' + window_name + '_EIVD_zeds.npy'
	np.save(new_EIVD_zed_file, EIVD_zeds)


	# Only column 6 is trait specific; every other path points at the shared artifact.
	t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + str(gwas_sample_size) + '\t')
	t.write(reference_rsid_file + '\t' + new_EIVD_zed_file + '\t' + EIVD_Q_file + '\t' + linked_genes_file + '\t' + snp_gene_anno_file + '\n')
f.close()
t.close()


###########
# Report
###########
print('##############################################################')
print('windows used: ' + str(n_windows_used) + '  |  skipped: ' + str(n_windows_skipped))
print('regression snps: ' + str(n_snps_total) + '  |  imputed: ' + str(n_snps_imputed))
if len(imputation_r2_arr) > 0:
	imputation_r2_arr = np.hstack(imputation_r2_arr)
	print('imputation r2: median ' + str(np.round(np.median(imputation_r2_arr),3)) + '  |  fraction below 0.5: ' + str(np.round(np.mean(imputation_r2_arr < 0.5),3)))
print('##############################################################')
