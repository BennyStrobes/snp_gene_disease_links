import numpy as np
import os
import sys
import pdb







def average_snp_gene_priors(full_output_stem, burn_in_iter, max_iter):
	prior_probs = []
	for itera in np.arange(burn_in_iter, max_iter+1, 5):
		file_name = full_output_stem + '_snp_gene_priors_iteration_' + str(itera) + '.txt'
		if os.path.exists(file_name) == False:
			continue
		data = np.loadtxt(file_name, dtype=str,delimiter='\t')
		prob_names = data[1:,0]
		probs = data[1:,1].astype(float)
		prior_probs.append(probs)
	prior_probs = np.asarray(prior_probs)

	output_file = full_output_stem + '_snp_gene_priors_averaged.txt'
	t = open(output_file,'w')
	t.write('prior_name\tprior_mean\tprior_mean_lb\tprior_mean_ub\n')
	for kk in range(prior_probs.shape[1]):

		meaner = np.mean(prior_probs[:,kk])

		lower_bound, upper_bound = np.percentile(prior_probs[:,kk], [2.5, 97.5])

		t.write(prob_names[kk] + '\t' + str(meaner) + '\t' + str(lower_bound) + '\t' + str(upper_bound) + '\n')

	t.close()

	return output_file



def average_gene_scores(full_output_stem, burn_in_iter, max_iter):
	gene_scores = []
	for ii,itera in enumerate(np.arange(burn_in_iter, max_iter+1, 5)):
		file_name = full_output_stem + '_gene_score_iteration_' + str(itera) + '.txt'
		if os.path.exists(file_name) == False:
			continue
		data = np.loadtxt(file_name, dtype=str,delimiter='\t')
		tmp_gene_names = data[1:,0]
		tmp_gene_scores = data[1:,1].astype(float)
		gene_scores.append(tmp_gene_scores)
		if ii == 0:
			gene_names = np.copy(tmp_gene_names)
		else:
			if np.array_equal(gene_names,tmp_gene_names) == False:
				print('assumption erororo')
				pdb.set_trace()

	gene_scores = np.asarray(gene_scores)

	output_file = full_output_stem + '_gene_score_averaged.txt'
	t = open(output_file,'w')
	t.write('gene_name\tgene_score\tgene_score_lb\tgene_score_ub\tgene_score_variance\tgene_z_score\n')

	means = []
	for gg, gene_name in enumerate(gene_names):
		vec = gene_scores[:, gg]

		if len(np.unique(vec)) == 1:
			continue

		meaner = np.mean(vec)
		means.append(meaner)
		lower_bound, upper_bound = np.percentile(vec, [2.5, 97.5])
		variance = np.var(vec)

		zs = meaner/np.sqrt(variance)



		t.write(gene_name + '\t' + str(meaner) + '\t' + str(lower_bound) + '\t' + str(upper_bound) + '\t' + str(variance) +'\t' + str(meaner/np.sqrt(variance)) + '\n')

	t.close()
	return output_file


trait_name = sys.argv[1]
output_stem = sys.argv[2]
prior_choice = sys.argv[3]
method_version = sys.argv[4]

# Relevent fields
full_output_stem = output_stem + '_lmm_snp_gene_link_' + prior_choice + '_' + method_version


burn_in_iter = 1500
max_iter = 4500


# Get average acrosss gibbs samples

# snp-gene priors
snp_gene_prior_output_file = average_snp_gene_priors(full_output_stem, burn_in_iter, max_iter)

# Gene scores output
gene_scores_output_file = average_gene_scores(full_output_stem, burn_in_iter, max_iter)

#print(gene_scores_output_file)

