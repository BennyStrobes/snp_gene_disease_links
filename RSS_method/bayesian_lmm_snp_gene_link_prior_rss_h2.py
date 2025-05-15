import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time
import scipy.special
from numba import njit

@njit(cache=True)
def fast_log_sum_exp_vector(aa):
	a_max = np.max(aa)                     # scalar
	shifted = aa - a_max                   # avoid overflow
	sum_exp = np.sum(np.exp(shifted))    # scalar
	return a_max + np.log(sum_exp)

def print_snp_gene_priors_output_file(pis, snp_gene_priors_output_file):
	t = open(snp_gene_priors_output_file,'w')
	t.write('component_name\tpi\n')
	for index, pi in enumerate(pis):
		t.write('component_' + str(index) + '\t' + str(pi) + '\n')

	t.close()
	return

@njit(cache=True)
def fast_categorical_sample(probs):
	cum_probs = np.cumsum(probs)
	rr = np.random.rand()
	return np.searchsorted(cum_probs, rr)


def print_gene_scores_to_output_file(gamma_vars, ordered_genes, gene_score_output_file):
	t = open(gene_score_output_file,'w')
	t.write('gene_name\tgene_score\n')
	for g_index, gene_name in enumerate(ordered_genes):
		t.write(gene_name + '\t' + str(gamma_vars[g_index]) + '\n')
	t.close()
	return

def print_snp_gene_link_output_file(window_names, class_membership, window_info, ordered_genes, snp_gene_link_output_file, method_version):
	t = open(snp_gene_link_output_file,'w')
	t.write('snp_id\tcategorical_class_membership\tnamed_class_membership\tgene_integer_class_membership\n')

	for window_name in window_names:
		window_class_membership = class_membership[window_name]
		window_snp_gene_mat = np.load(window_info[window_name]['snp_gene_names_file'])
		window_rsids = window_info[window_name]['snp_names']

		for snp_index, snp_name in enumerate(window_rsids):
			snp_class = window_class_membership[snp_index]
			if method_version.startswith('null_component'):
				if snp_class == 0:
					t.write(snp_name + '\t' + str(snp_class) + '\t' 'NULL' + '\t' + 'NULL' + '\n')
				else:
					gene_index = window_snp_gene_mat[snp_index, (snp_class-1)]
					gene_name = ordered_genes[gene_index]
					t.write(snp_name + '\t' + str(snp_class) + '\t' + gene_name + '\t' + str(gene_index) + '\n')
			elif method_version.startswith('snp_gene_component'):
				gene_index = window_snp_gene_mat[snp_index, (snp_class)]
				gene_name = ordered_genes[gene_index]
				t.write(snp_name + '\t' + str(snp_class) + '\t' + gene_name + '\t' + str(gene_index) + '\n')			
	t.close()
	return


def update_gamma_from_single_window(QQ, gamma_vec, gwas_beta_resid, gamma_var, N_gwas, resid_var, pis, null_variance, method_version, snp_gene_link_mat, univariate=True):

	# N-snps in window
	KK = len(gamma_vec)


	gamma_variance_grid = gamma_var[snp_gene_link_mat]

	if method_version == 'null_component':
		gamma_variance_grid = np.hstack((null_variance*np.ones((gamma_variance_grid.shape[0], 1)), gamma_variance_grid))


	# Precompute posterior variance (same for all snps)
	t0 = time.time()
	qq_sq = np.sum(np.square(QQ),axis=0)

	class_membership_vec = np.zeros(KK).astype(int)

	posterior_var_term1 = (qq_sq*N_gwas/resid_var)
	gamma_precision_grid = 1.0/gamma_variance_grid

	#t15 = time.time()
	posterior_var_grid = 1.0 / (posterior_var_term1[:, np.newaxis] + gamma_precision_grid)
	#t16 = time.time()


	# Update each snp in turn
	for snp_index in np.random.permutation(range(KK)):

		# Range of posterior vars
		#posterior_var = 1.0/((qq_sq[snp_index]*N_gwas/resid_var) + (1.0/gamma_variance_grid[snp_index, :]))
		#posterior_var = 1.0/(posterior_var_term1[snp_index] + gamma_precision_grid[snp_index,:])
		posterior_var = posterior_var_grid[snp_index,:]

		# Re include current effect
		gwas_beta_resid = gwas_beta_resid + (QQ[:, snp_index]*gamma_vec[snp_index])

		# Compute posterior mean for this snp
		posterior_mean = posterior_var*(N_gwas/resid_var)*np.dot(gwas_beta_resid,QQ[:,snp_index])

		# Compute log likelihoods of each class classes
		log_like = np.log(pis) - (.5*np.log(gamma_variance_grid[snp_index,:]/posterior_var)) + (.5*(np.square(posterior_mean)/(posterior_var)))

		# Fix zeros
		zero_components = posterior_var == 0.0
		log_like[zero_components] = np.log(pis[zero_components])

		# Get probability of each class
		temp = -log_like + fast_log_sum_exp_vector(log_like)
		#temp = -log_like + scipy.special.logsumexp(log_like)
		#temp = scipy.special.logsumexp(log_like - log_like[:, np.newaxis],axis=1)
		#if np.sum(temp >600) > 0:
			#temp[temp>600] = 600
		probs = 1.0/np.exp(temp)

		# Sample class membership
		class_p = fast_categorical_sample(probs)
		#t8 = time.time()
		#class_p = np.random.choice(np.arange(len(probs)), p=probs)
		class_membership_vec[snp_index] = class_p

		# Sample from posterior distribution
		gamma_vec[snp_index] = np.random.normal(loc=posterior_mean[class_p], scale=np.sqrt(posterior_var[class_p]))

		# Remove updated effect
		gwas_beta_resid = gwas_beta_resid - (QQ[:, snp_index]*gamma_vec[snp_index])

	return gamma_vec, gwas_beta_resid, class_membership_vec





class Bayesian_LMM_RSS_h2_inference(object):
	def __init__(self, window_names, window_info, n_gwas_individuals, ordered_genes, output_stem,inv_gamma_alpha=1e-10, inv_gamma_beta=1e-10, cross_gene_hyperparm=False,method_version='null_component'):
		self.N_gwas = n_gwas_individuals
		self.window_names = window_names
		self.window_info = window_info
		self.ordered_genes = np.asarray(ordered_genes)
		self.output_stem = output_stem
		self.cross_gene_hyperparm = cross_gene_hyperparm

		self.inv_gamma_alpha = inv_gamma_alpha
		self.inv_gamma_beta = inv_gamma_beta
		self.method_version = method_version

		# Number of snps
		self.KK = 0.0
		for window in window_names:
			self.KK = self.KK + self.window_info[window]['n_snps']
		
		return



	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_resid_var_bool=True):
		""" Fit the model.
		"""
		# Initialize model params
		self.initialize_variables()

		# Keep track of iterations
		self.itera = 0

		#print(self.output_log_file)
		t = open(self.output_stem + 'log.txt','w')
		print(self.output_stem + 'log.txt')

		# Iterative Gibbs sampling algorithm
		for itera in range(total_iterations):
			print(itera)
			
			# Update gamma
			t1 = time.time()
			self.update_gamma()
			t2 = time.time()
			print(t2-t1)

			# Update PI
			if self.method_version == 'snp_gene_component_smart_init':
				if itera > 5:
					self.update_pi()
			elif self.method_version == 'snp_gene_component_fixed_to_smart_init':
				t1 = time.time()
			else:
				self.update_pi()
	
			# Update gamma_var
			self.update_gamma_var()

			if self.cross_gene_hyperparm:
				self.update_inv_gamma_beta_prior()

			# Update resid var
			if update_resid_var_bool:
				self.update_resid_var()



			# Update iteration number
			self.itera = self.itera + 1

			t.write('Iteration: ' + str(itera) + '\n')
			t.write(','.join(self.pis.astype(str)) + '\n')
			t.write(str(self.resid_var) + '\n')
			t.write(str(self.inv_gamma_beta/self.inv_gamma_alpha) + '\n')
			t.flush()
			print(','.join(self.pis.astype(str)))


			if np.mod(self.itera, 5) == 0:
				# Gene scores ouput
				gene_score_output_file = self.output_stem + 'gene_score_iteration_' + str(self.itera) + '.txt'
				print_gene_scores_to_output_file(self.gamma_var, self.ordered_genes, gene_score_output_file)
				# Snp gene link output
				snp_gene_link_output_file = self.output_stem + 'snp_gene_links_iteration_' + str(self.itera) + '.txt'
				print_snp_gene_link_output_file(self.window_names, self.class_membership, self.window_info, self.ordered_genes, snp_gene_link_output_file, self.method_version)
				# Snp gene priors output
				snp_gene_priors_output_file = self.output_stem + 'snp_gene_priors_iteration_' + str(self.itera) + '.txt'
				print_snp_gene_priors_output_file(self.pis, snp_gene_priors_output_file)

			'''
			if np.mod(self.itera,100) == 0.0:
				print(self.gamma_var*self.KK)
			'''
			#t.write('Iteration: ' + str(itera) + '\n')
			#t.write(','.join(self.pis.astype(str)) + '\n')
			#t.write(str(np.sum(self.gamma_var*self.deltas*self.pis*self.KK)) + '\n')
			#t.write(str(self.resid_var) + '\n')
			#t.flush()

			if itera > burn_in_iterations:
				self.sampled_h2.append(self.gamma_var*self.KK)
				self.sampled_resid_var.append(self.resid_var)


		self.sampled_h2 = np.asarray(self.sampled_h2)
		self.sampled_resid_var = np.asarray(self.sampled_resid_var)
		t.close()

		return

	def update_pi(self):
		counts = np.zeros(self.C)
		for window_name in self.window_names:
			window_membership_vec = self.class_membership[window_name]

			#window_snp_gene_mat = np.load(self.window_info[window_name]['snp_gene_names_file'])

			#linked_genes = window_snp_gene_mat[np.arange(window_snp_gene_mat.shape[0]),window_membership_vec]

			#valid_indices = self.gamma_var[linked_genes] > 1e-10

			for class_name in np.arange(self.C):
				counts[class_name] = counts[class_name] + np.sum(window_membership_vec ==class_name)

		# Randomly sample pi from dirichlet
		self.pis = np.random.dirichlet(counts + self.alpha_0)

		# Set values of zero to really small number
		self.pis[self.pis==0.0]=1e-30

		return


	def update_resid_var(self, v0=0.0, s_sq=0.0, cc=1e-8):
		all_resids = []
		for window_name in self.window_names:
			window_resid = self.beta_resid[window_name]
			all_resids.append(window_resid)
		all_resids = np.hstack(all_resids)

		vv = len(all_resids) + v0
		tau_sq = np.sum(np.square(all_resids)*self.N_gwas) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		self.resid_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_inv_gamma_beta_prior(self, beta0_tmp=0.0, alpha0_tmp=0.0):
		
		gamma_b =beta0_tmp + np.sum(1.0/self.gamma_var[self.valid_gene_indices])
		gamma_a = alpha0_tmp + len(self.gamma_var[self.valid_gene_indices])*self.inv_gamma_alpha

		self.inv_gamma_beta = np.random.gamma(gamma_a, scale=1/gamma_b,size=1)[0]

		return

	def update_gamma_var(self, v0=0.0, s_sq=0.0):
		# First get middle gammas
		sum_gamma_sq_vec = np.zeros(self.G)
		num_snps = np.zeros(self.G)
		for window_name in self.window_names:
			window_sq_gammas = np.square(self.gamma[window_name])
			window_class_membership = self.class_membership[window_name]
			window_snp_gene_mat = np.load(self.window_info[window_name]['snp_gene_names_file'])

			if self.method_version == 'null_component':
				gene_linked_indices = window_class_membership != 0.0
				gene_linked_window_class_membership = window_class_membership[gene_linked_indices] - 1
				gene_linked_window_snp_gene_mat = window_snp_gene_mat[gene_linked_indices, :]
				linked_genes = gene_linked_window_snp_gene_mat[np.arange(gene_linked_window_snp_gene_mat.shape[0]),gene_linked_window_class_membership]
				np.add.at(sum_gamma_sq_vec, linked_genes, window_sq_gammas[gene_linked_indices])
				np.add.at(num_snps, linked_genes, np.ones(len(window_sq_gammas[gene_linked_indices])))
			elif self.method_version == 'snp_gene_component' or self.method_version == 'snp_gene_component_smart_init' or self.method_version == 'snp_gene_component_fixed_to_smart_init':
				linked_genes = window_snp_gene_mat[np.arange(window_snp_gene_mat.shape[0]),window_class_membership]
				np.add.at(sum_gamma_sq_vec, linked_genes, window_sq_gammas)
				np.add.at(num_snps, linked_genes, np.ones(len(window_sq_gammas)))




		vv = num_snps + v0
		tau_sq = sum_gamma_sq_vec + s_sq

		self.valid_gene_indices = num_snps >= 1


		param1 = (vv/2) + self.inv_gamma_alpha # inv_gamma_alpha == gene_gamma_shape
		param2 = 1.0/((tau_sq/2) + self.inv_gamma_beta) # inv_gamma_beta == self.gene_gamma_scale)


		self.gamma_var[self.valid_gene_indices] = 1.0/np.random.gamma(shape=param1[self.valid_gene_indices], scale=param2[self.valid_gene_indices])

		return

	def update_gamma(self):
		for window_name in self.window_names:
			gamma_vec, gwas_beta_resid_vec, class_membership_vec = update_gamma_from_single_window(np.load(self.window_info[window_name]['Q_file']), self.gamma[window_name], self.beta_resid[window_name], self.gamma_var, self.N_gwas, self.resid_var, self.pis, self.null_variance, self.method_version, np.load(self.window_info[window_name]['snp_gene_names_file']))
			self.gamma[window_name] = gamma_vec
			self.beta_resid[window_name] = gwas_beta_resid_vec
			self.class_membership[window_name] = class_membership_vec

		return



	def initialize_variables(self):
		# Initialize causal effcts
		self.gamma = {}
		self.beta_resid = {}
		self.class_membership = {}

		for window_name in self.window_names:
			self.gamma[window_name] = np.zeros(self.window_info[window_name]['n_snps'])
			self.beta_resid[window_name] = np.copy(self.window_info[window_name]['beta_pc'])
			self.class_membership[window_name] = np.zeros(self.window_info[window_name]['n_snps']).astype(int)


		random_snp_gene_names_mat = np.load(self.window_info[window_name]['snp_gene_names_file'])

		# Number of categories
		if self.method_version == 'null_component':
			self.C = random_snp_gene_names_mat.shape[1] + 1
		elif self.method_version == 'snp_gene_component' or self.method_version == 'snp_gene_component_smart_init' or self.method_version == 'snp_gene_component_fixed_to_smart_init':
			self.C = random_snp_gene_names_mat.shape[1]

		self.null_variance=1e-20

		# Number of genes
		self.G = len(self.ordered_genes)

		# Initialize variance parameters
		self.gamma_var = np.ones(self.G)*.1/1000000

		# Initialize residual variance
		self.resid_var = 1.0

		# Initialize pis
		self.pis = np.ones(self.C)/self.C
		if self.method_version == 'snp_gene_component_smart_init' or self.method_version == 'snp_gene_component_fixed_to_smart_init':
			self.pis[0] = .314
			self.pis[1] = .1279
			self.pis[2] = .0879
			self.pis[3] = .0879
			self.pis[4] = .0879
			self.pis[5] = .0098
			self.pis[6] = .0098
			self.pis[7] = .0098
			self.pis[8] = .0098
			self.pis[9] = .0098
			self.pis = self.pis/np.sum(self.pis)


		# Initialize hyperparameter on pis
		self.alpha_0 = np.ones(self.C)

		# Keep track of sampled gamma_vars
		self.sampled_h2 = []
		self.sampled_resid_var = []

		return


