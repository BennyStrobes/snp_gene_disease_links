import sys
import numpy as np 
import os
import pdb
import bayesian_lmm_snp_gene_link_prior_rss_h2

def count_rsids_in_window(rsid_file):
	"""Number of snps in a window's rsid file.

	Only the COUNT is needed while the chain runs (it sizes gamma). The rsid strings
	themselves are re-read once at output time -- holding all ~9.9M of them for the whole
	run costs ~435 MB for something used only when writing the final links file.
	"""
	n_snps = 0
	f = open(rsid_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			head_count = head_count + 1
			continue
		n_snps = n_snps + 1
	f.close()

	return n_snps

def load_in_data(input_window_summary_file):
	window_names = []
	window_Q_files = []
	window_rsid_files = []
	window_n_snps = []
	window_zs = []
	window_snp_gene_names = []
	window_snp_gene_annotations = []

	gwas_sample_size = None

	head_count = 0
	f = open(input_window_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields
		window_name = data[0]
		rsid_file = data[5]
		zscore_file = data[6]
		Q_mat_file = data[7]
		linked_genes_file = data[8]
		snp_gene_anno_file = data[9]

		# The gwas sample size is a property of the trait, not of the window
		line_gwas_sample_size = float(data[4])
		if gwas_sample_size is None:
			gwas_sample_size = line_gwas_sample_size
		elif line_gwas_sample_size != gwas_sample_size:
			print('assumption eroror: gwas sample size varies across windows (' + str(gwas_sample_size) + ' vs ' + str(line_gwas_sample_size) + ')')
			pdb.set_trace()

		n_snps = count_rsids_in_window(rsid_file)
		zs = np.load(zscore_file)

		# Under asymmetric LD Q_mat is (n_whitened_components X n_reference_snps) and is NOT
		# square. Its COLUMNS are the reference snps -- they carry gamma and the gene links.
		# Its ROWS are the whitened regression-snp space -- they carry the gwas z-scores.
		# Both axes are checked here because a half-migrated pipeline otherwise fails deep
		# inside the gibbs sampler with an opaque index error.
		# mmap_mode reads only the .npy header, so this costs no memory.
		Q_shape = np.load(Q_mat_file, mmap_mode='r').shape
		linked_genes_shape = np.load(linked_genes_file, mmap_mode='r').shape

		if Q_shape[1] != n_snps:
			print('assumption eroror: window ' + window_name + ' Q_mat has ' + str(Q_shape[1]) + ' columns but ' + str(n_snps) + ' snps in ' + rsid_file)
			print('Column 5 of the window summary must be the REFERENCE snp file, not the regression snp file.')
			pdb.set_trace()
		if linked_genes_shape[0] != n_snps:
			print('assumption eroror: window ' + window_name + ' linked genes has ' + str(linked_genes_shape[0]) + ' rows but ' + str(n_snps) + ' reference snps')
			pdb.set_trace()
		if Q_shape[0] != len(zs):
			print('assumption eroror: window ' + window_name + ' Q_mat has ' + str(Q_shape[0]) + ' rows but ' + str(len(zs)) + ' whitened z-scores')
			pdb.set_trace()

		# Add to global arrays
		window_names.append(window_name)
		window_rsid_files.append(rsid_file)
		window_n_snps.append(n_snps)
		window_Q_files.append(Q_mat_file)
		window_zs.append(zs)
		window_snp_gene_names.append(linked_genes_file)
		window_snp_gene_annotations.append(snp_gene_anno_file)

	f.close()

	if gwas_sample_size is None:
		print('assumption eroror: no windows found in ' + input_window_summary_file)
		pdb.set_trace()

	print('windows: ' + str(len(window_names)) + ' | gwas sample size: ' + str(gwas_sample_size))
	print('reference snps (gamma): ' + str(np.sum(window_n_snps)) + ' | whitened components (z): ' + str(np.sum([len(x) for x in window_zs])))

	return np.asarray(window_names), window_rsid_files, window_n_snps, window_zs, np.asarray(window_Q_files), np.asarray(window_snp_gene_names), np.asarray(window_snp_gene_annotations), gwas_sample_size

def load_in_genes(gene_summary_file):
	gene_arr = []
	f = open(gene_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_arr.append(data[0])
	f.close()
	return np.asarray(gene_arr)


def write_gene_scores_to_output(gene_score_output_file, gene_scores, ordered_genes):
	t = open(gene_score_output_file,'w')
	t.write('gene_name\tgene_score\n')
	for gene_iter, gene_name in enumerate(ordered_genes):
		t.write(gene_name + '\t' + str(gene_scores[gene_iter]) + '\n')
	t.close()
	return

def write_snp_gene_weight_parameters_to_output(snp_gene_weight_parameters_output_file, intercept, annotation_weights):
	t = open(snp_gene_weight_parameters_output_file,'w')
	t.write('parameter_name\tweight\n')
	t.write('intercept\t' + str(intercept) + '\n')

	for ii, weight in enumerate(annotation_weights):
		t.write('anno_' + str(ii) + '\t' + str(weight) + '\n')
	t.close()
	return

def write_average_gene_rank_probabilities(average_gene_rank_probabilities_output_file, gene_rank_probs):
	new_probs = gene_rank_probs/np.sum(gene_rank_probs)
	t = open(average_gene_rank_probabilities_output_file,'w')
	t.write('gene_rank\tprobability\n')
	for gene_iter, gene_probs in enumerate(new_probs):
		t.write('rank_' + str(gene_iter+1) + '\t' + str(gene_probs) + '\n')
	t.close()
	return

#######################
# Command line args
#######################
trait_name = sys.argv[1]
input_window_summary_file = sys.argv[2]
gene_summary_file = sys.argv[3]
output_stem = sys.argv[4]
prior_choice = sys.argv[5]
method_version = sys.argv[6]

################
# Learning parameters
# Only one method version is supported; everything else has been removed
if method_version != 'snp_gene_component_fixed_to_smart_init':
	print('assumption eroror: only snp_gene_component_fixed_to_smart_init is supported, got ' + str(method_version))
	pdb.set_trace()

######
# Load in genes
ordered_genes = load_in_genes(gene_summary_file)



######
# Load in data
# 1. Window names
# 2. Window rsids
# 3. Window z-scores (in PC space)
# 4. Window Q files (basically LD equivalent in PC space)
# 5. Window snp gene names 
# 6. Window snp gene annos 
# 7. GWAS sample size
window_names, window_rsid_files, window_n_snps, window_zs, window_Q_files, window_snp_gene_name_files, window_snp_gene_annotation_files, gwas_sample_size = load_in_data(input_window_summary_file)

##*#*#*#*#*#**#
#window_names = window_names[:130]
##*#*#*#*#*#**#
#num_windows = len(window_names)


#############
# First Run standard bayesian heritability model
# Organize the data
window_info = {}
for ii, window_name in enumerate(window_names):
	window_info[window_name] = {}
	window_info[window_name]['n_snps'] = window_n_snps[ii]
	window_info[window_name]['rsid_file'] = window_rsid_files[ii]
	window_info[window_name]['Q_file'] = window_Q_files[ii]
	window_info[window_name]['beta_pc'] = window_zs[ii]/np.sqrt(gwas_sample_size)
	window_info[window_name]['snp_gene_names_file'] = window_snp_gene_name_files[ii]



if prior_choice.startswith('inverse_gamma_1e-'):
	inv_gamma_alpha_prior = float(prior_choice.split('_gamma_')[1])
	inv_gamma_beta_prior = float(prior_choice.split('_gamma_')[1])

if prior_choice.startswith('inverse_gamma_cross_gene_prior'):
	cross_gene_hyperparm=True
	inv_gamma_alpha_prior = float(prior_choice.split('cross_gene_prior_')[1])
	initial_value = 1e-7
	inv_gamma_beta_prior = initial_value*inv_gamma_alpha_prior
else:
	cross_gene_hyperparm = False


mod = bayesian_lmm_snp_gene_link_prior_rss_h2.Bayesian_LMM_RSS_h2_inference(window_names, window_info, gwas_sample_size, ordered_genes, output_stem + '_lmm_snp_gene_link_' + prior_choice + '_' + method_version + '_', inv_gamma_alpha=inv_gamma_alpha_prior, inv_gamma_beta=inv_gamma_beta_prior, method_version=method_version, cross_gene_hyperparm=cross_gene_hyperparm)
mod.fit(total_iterations=7000, burn_in_iterations=5000, update_resid_var_bool=True, thin_iterations=2)

















'''
mod = bayesian_lmm_rss_h2.Bayesian_LMM_RSS_h2_inference(window_names, window_info, gwas_sample_size, output_stem + '_lmm_log.txt')
mod.fit()
'''

'''
mod = bayesian_lmm_mog_prior_rss_h2.Bayesian_LMM_RSS_h2_inference(window_names, window_info, gwas_sample_size, output_stem + '_lmm_mog_log.txt')
mod.fit()
'''



'''
#OLD!
# Extract number of features
tmp_window_snp_gene_anno = np.load(window_snp_gene_annotation_files[0])
n_samp, n_nearby_genes, K = tmp_window_snp_gene_anno.shape


optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)

W = tf.Variable(1e-8*tf.random.normal(shape=(K,)), dtype=tf.float32) # Weights
bb = tf.Variable(0.0, dtype=tf.float32)  # (K+1,)
cc = tf.Variable(0.0, dtype=tf.float32)  # (K+1,)
sig_sq_g = tf.Variable(1e-16*tf.abs(tf.random.normal(shape=(len(ordered_genes),))), dtype=tf.float32) # Weights


for epoch_iter in range(max_epochs):
	print('###################################')
	print('epoch iter ' + str(epoch_iter))
	print('###################################')
	if null_component:
		gene_rank_probs = np.zeros(n_nearby_genes+1)
	else:
		gene_rank_probs = np.zeros(n_nearby_genes)
	# Loop through genomic windows
	#for window_counter, window_iter in enumerate(np.random.permutation(range(num_windows))):
	for window_counter, window_iter in enumerate(range(num_windows)):

		# Get data for this window
		window_name = window_names[window_iter]


		window_sq_ld_file = window_sq_ld_files[window_iter]
		window_regression_snp_filter = window_regression_snp_filters[window_iter]
		window_chi_sq = window_chi_sqs[window_iter]
		window_snp_gene_name_file = window_snp_gene_name_files[window_iter]
		window_snp_gene_annotation_file = window_snp_gene_annotation_files[window_iter]

		# Load in data
		window_sq_ld = np.load(window_sq_ld_file)[window_regression_snp_filter,:]
		window_snp_gene_name = np.load(window_snp_gene_name_file)
		window_snp_gene_anno = np.load(window_snp_gene_annotation_file)

		# Dimensionality of space
		n_samp, n_nearby_genes, K = window_snp_gene_anno.shape

		with tf.GradientTape() as tape:
			# predict per snp heritabilities
			per_snp_h2, window_gene_rank_probs = pred_per_snp_h2(W, bb, cc, sig_sq_g, window_snp_gene_anno, window_snp_gene_name, n_samp, scaling_factor, epoch_iter, null_component=null_component)

			# Compute Loss
			loss_value = ldsc_tf_loss_fxn(window_chi_sq, per_snp_h2, gwas_sample_size, window_sq_ld)

		#pdb.set_trace()
		#print('loss')
		#print(loss_value)

		# Compute gradients
		#gradients = tape.gradient(loss_value, [W, bb, cc, sig_sq_g])
		gradients = tape.gradient(loss_value, [W, sig_sq_g])
		
		# Update weights using gradient descent
		#optimizer.apply_gradients(zip(gradients, [W, bb, cc, sig_sq_g]))
		optimizer.apply_gradients(zip(gradients, [W, sig_sq_g]))

		gene_rank_probs = gene_rank_probs + window_gene_rank_probs

		if np.mod(window_counter,100) == 0.0:
			print(window_counter)
			print(W)
			print(np.sort(tf.math.exp(sig_sq_g)/scaling_factor))
		#print(np.sort(sig_sq_g/scaling_factor))

	# Print files to output
	gene_score_output_file = output_stem + '_gene_scores_' + str(epoch_iter) + '.txt'
	write_gene_scores_to_output(gene_score_output_file, tf.math.exp(sig_sq_g).numpy()/scaling_factor, ordered_genes)

	# Print files to output
	snp_gene_weight_parameters_output_file = output_stem + '_snp_gene_weight_parameters_' + str(epoch_iter) + '.txt'
	write_snp_gene_weight_parameters_to_output(snp_gene_weight_parameters_output_file, bb.numpy(), W.numpy())

	# Print gene rank probabilities to output file
	average_gene_rank_probabilities_output_file = output_stem + '_average_gene_rank_probabilities_' + str(epoch_iter) + '.txt'
	write_average_gene_rank_probabilities(average_gene_rank_probabilities_output_file, gene_rank_probs)
	print(average_gene_rank_probabilities_output_file)
'''


# logits = tf.squeeze(tf.matmul(np.float32(window_snp_gene_anno), W) +bb,axis=-1)
# new_column = tf.fill([n_samp, 1], cc)
# logits2 = tf.concat([new_column, logits], axis=1)
# probs = tf.nn.softmax(logits2,axis=1)
