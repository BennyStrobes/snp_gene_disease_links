import sys
import numpy as np 
import pandas as pd
import os
import pdb
import tensorflow as tf
import gzip
import time
from tensorflow.keras.layers import Conv1D, Dense, BatchNormalization, Activation, GlobalAveragePooling1D, Add, Input, Flatten
from tensorflow.keras import Model





def create_dictionary_mapping_from_rsid_to_chi_sq(trait_sumstat_file):
	f = open(trait_sumstat_file)
	sample_size = None
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[0]
		if sample_size is None:
			sample_size = int(data[3])
		else:
			if sample_size != int(data[3]):
				print('assumption eroror')
				pdb.set_trace()

		chi_sq = float(data[4])*float(data[4])

		if rsid in dicti:
			print('assumption erororo')
			pdb.set_trace()
		dicti[rsid] = chi_sq

	f.close()


	return dicti, sample_size


def extract_middle_regression_rsids(middle_regression_rsid_file):
	ll = open(middle_regression_rsid_file)
	rsids = []
	head_count = 0
	for line in ll:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsids.append(line)
	ll.close()

	return np.asarray(rsids)

def get_chi_sq_stats_for_this_window(middle_regression_rsids, rsid_to_chi_sq):
	chi_sq_arr = []
	filter_arr = []

	for rsid in middle_regression_rsids:
		if rsid in rsid_to_chi_sq:
			chi_sq_arr.append(rsid_to_chi_sq[rsid])
			filter_arr.append(True)
		else:
			filter_arr.append(False)
			#print('miss')


	chi_sq_arr = np.asarray(chi_sq_arr)
	filter_arr = np.asarray(filter_arr)

	if len(chi_sq_arr) < 2:
		print('assumption eroror')
		pdb.set_trace()

	return chi_sq_arr, filter_arr


def load_in_data(input_window_summary_file, rsid_to_chi_sq):
	window_names = []
	window_ld_files = []
	window_regression_snp_filters = []
	window_chi_sqs = []
	window_snp_gene_names = []
	window_snp_gene_annotations = []

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
		LD_file = data[4]
		middle_regression_rsid_file = data[6]
		linked_genes_file = data[7]
		snp_gene_anno_file = data[8]

		# Extract middle regression rsids
		middle_regression_rsids = extract_middle_regression_rsids(middle_regression_rsid_file)

		# Get chi-sq stats for this window
		window_chi_sq, window_regression_snp_filter = get_chi_sq_stats_for_this_window(middle_regression_rsids, rsid_to_chi_sq)
		
		# Get linkded genes
		#linked_genes = np.load(linked_genes_file)

		# Add to global arrays
		window_names.append(window_name)
		window_ld_files.append(LD_file)
		window_regression_snp_filters.append(window_regression_snp_filter)
		window_chi_sqs.append(window_chi_sq)
		window_snp_gene_names.append(linked_genes_file)
		window_snp_gene_annotations.append(snp_gene_anno_file)


	f.close()



	return np.asarray(window_names), np.asarray(window_ld_files), window_regression_snp_filters, window_chi_sqs, np.asarray(window_snp_gene_names), np.asarray(window_snp_gene_annotations)

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

def pred_per_snp_h2(W, bb, cc, sig_sq_g, window_snp_gene_anno, window_snp_gene_name, n_samp):
	#logits = tf.squeeze(tf.matmul(np.float32(window_snp_gene_anno), W) +bb,axis=-1)
	logits = tf.linalg.matvec(np.float32(window_snp_gene_anno), W) + bb

	if np.sum(tf.math.is_nan(logits)) > 0:
		pdb.set_trace()

	#new_column = tf.fill([n_samp, 1], cc)
	#logits2 = tf.concat([new_column, logits], axis=1)
	#probs = tf.nn.softmax(logits2,axis=1)
	probs = tf.nn.softmax(logits,axis=1)

	if np.sum(tf.math.is_nan(probs)) > 0:
		pdb.set_trace()

	matrix_sig_sq_g = tf.gather(sig_sq_g, window_snp_gene_name)

	if np.sum(tf.math.is_nan(matrix_sig_sq_g)) > 0:
		pdb.set_trace()

	#new_column2 = tf.fill([n_samp, 1], 0.0)
	#matrix_sig_sq_g2 = tf.concat([new_column2, matrix_sig_sq_g], axis=1)

	#per_snp_h2 = tf.reduce_sum(tf.multiply(probs, matrix_sig_sq_g2), axis=1)
	per_snp_h2 = tf.reduce_sum(tf.multiply(probs, matrix_sig_sq_g/100000.0), axis=1)

	print('per snp h2')
	print(np.sort(np.abs(per_snp_h2)))


	return per_snp_h2

def ldsc_tf_loss_fxn(window_chi_sq, per_snp_h2, gwas_sample_size, window_sq_ld):
	pred_chi_sq = (gwas_sample_size*tf.linalg.matvec(np.float32(window_sq_ld), per_snp_h2)) + 1.0
	squared_error = tf.pow(pred_chi_sq - window_chi_sq,2)

	if np.sum(tf.math.is_nan(squared_error)) > 0:
		pdb.set_trace()

	if tf.math.reduce_sum(squared_error) > 1e15:
		pdb.set_trace()

	return tf.math.reduce_sum(squared_error)



#######################
# Command line args
#######################
trait_name = sys.argv[1]
trait_sumstat_file = sys.argv[2]
input_window_summary_file = sys.argv[3]
gene_summary_file = sys.argv[4]
output_stem = sys.argv[5]

################
# Learning parameters
learning_rate=1e-8 
max_epochs = 200


######
# Load in genes
ordered_genes = load_in_genes(gene_summary_file)

######
# Load in sumstats
# Create dictionary mapping from rsid to chi-squared statistic
# Also extract gwas sample size
rsid_to_chi_sq, gwas_sample_size = create_dictionary_mapping_from_rsid_to_chi_sq(trait_sumstat_file)

######
# Load in data
# 1. Window names
# 2. Window LD files
# 3. Window chi-squared stats
# 4. Window snp gene names (either data or files depending on mem)
# 5. Window snp gene annos (either data or files depending on mem)
window_names, window_sq_ld_files, window_regression_snp_filters, window_chi_sqs, window_snp_gene_name_files, window_snp_gene_annotation_files = load_in_data(input_window_summary_file, rsid_to_chi_sq)
num_windows = len(window_names)

# Extract number of features
tmp_window_snp_gene_anno = np.load(window_snp_gene_annotation_files[0])
n_samp, n_nearby_genes, K = tmp_window_snp_gene_anno.shape


optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)

W = tf.Variable(1e-8*tf.random.normal(shape=(K,)), dtype=tf.float32) # Weights
bb = tf.Variable(0.0, dtype=tf.float32)  # (K+1,)
cc = tf.Variable(0.0, dtype=tf.float32)  # (K+1,)
sig_sq_g = tf.Variable(1e-16*tf.abs(tf.random.normal(shape=(len(ordered_genes),))), dtype=tf.float32) # Weights

print(W)
print(sig_sq_g)

for epoch_iter in range(max_epochs):
	print('###################################')
	print('epoch iter ' + str(epoch_iter))
	print('###################################')

	# Loop through genomic windows
	for window_counter, window_iter in enumerate(np.random.permutation(range(num_windows))):

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
			per_snp_h2 = pred_per_snp_h2(W, bb, cc, sig_sq_g, window_snp_gene_anno, window_snp_gene_name, n_samp)

			# Compute Loss
			loss_value = ldsc_tf_loss_fxn(window_chi_sq, per_snp_h2, gwas_sample_size, window_sq_ld)

		#pdb.set_trace()
		print('loss')
		print(loss_value)

		# Compute gradients
		gradients = tape.gradient(loss_value, [W, bb, cc, sig_sq_g])
		
		# Update weights using gradient descent
		optimizer.apply_gradients(zip(gradients, [W, bb, cc, sig_sq_g]))

		print(W)
		print(sig_sq_g)



# logits = tf.squeeze(tf.matmul(np.float32(window_snp_gene_anno), W) +bb,axis=-1)
# new_column = tf.fill([n_samp, 1], cc)
# logits2 = tf.concat([new_column, logits], axis=1)
# probs = tf.nn.softmax(logits2,axis=1)
