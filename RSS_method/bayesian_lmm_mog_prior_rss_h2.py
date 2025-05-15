import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time
import scipy.special

def update_gamma_from_single_window(QQ, gamma_vec, gwas_beta_resid, gamma_var, N_gwas, resid_var, pis, deltas, univariate=True):

	# N-snps in window
	KK = len(gamma_vec)

	gamma_variance_grid = deltas*gamma_var

	# Precompute posterior variance (same for all snps)
	qq_sq = np.sum(np.square(QQ),axis=0)

	class_membership_vec = np.zeros(KK).astype(int)


	# Update each snp in turn
	for snp_index in np.random.permutation(range(KK)):

		# Range of posterior vars
		posterior_var = 1.0/((qq_sq[snp_index]*N_gwas/resid_var) + (1.0/gamma_variance_grid))

		# Re include current effect
		gwas_beta_resid = gwas_beta_resid + (QQ[:, snp_index]*gamma_vec[snp_index])
		# Compute posterior mean for this snp
		posterior_mean = posterior_var*(N_gwas/resid_var)*np.dot(gwas_beta_resid,QQ[:,snp_index])

		# Compute log likelihoods of each class classes
		log_like = np.log(pis) - (.5*np.log(gamma_variance_grid/posterior_var)) + (.5*(np.square(posterior_mean)/(posterior_var)))

		# Fix zeros
		zero_components = posterior_var == 0.0
		log_like[zero_components] = np.log(pis[zero_components])

		# Get probability of each class
		temp = scipy.special.logsumexp(log_like - log_like[:, np.newaxis],axis=1)
		if np.sum(temp >600):
			temp[temp>600] = 600
		probs = 1.0/np.exp(temp)


		# Sample class membership
		class_p = np.random.choice(np.arange(len(probs)), p=probs)
		class_membership_vec[snp_index] = class_p

		# Sample from posterior distribution
		gamma_vec[snp_index] = np.random.normal(loc=posterior_mean[class_p], scale=np.sqrt(posterior_var[class_p]))

		# Remove updated effect
		gwas_beta_resid = gwas_beta_resid - (QQ[:, snp_index]*gamma_vec[snp_index])


	return gamma_vec, gwas_beta_resid, class_membership_vec





class Bayesian_LMM_RSS_h2_inference(object):
	def __init__(self, window_names, window_info, n_gwas_individuals, output_log_file, deltas=[1e-20, .01, .1, 1.0]):
		self.N_gwas = n_gwas_individuals
		self.window_names = window_names
		self.window_info = window_info
		self.deltas = np.asarray(deltas)
		self.output_log_file = output_log_file

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
		print(self.output_log_file)

		t = open(self.output_log_file,'w')

		# Iterative Gibbs sampling algorithm
		for itera in range(total_iterations):
			print(itera)
			# Update gamma
			self.update_gamma()

			# Update PI
			self.update_pi()
	
			# Update gamma_var
			self.update_gamma_var()

			# Update resid var
			if update_resid_var_bool:
				self.update_resid_var()

			# Update iteration number
			self.itera = self.itera + 1

			'''
			if np.mod(self.itera,100) == 0.0:
				print(self.gamma_var*self.KK)
			'''
			t.write('Iteration: ' + str(itera) + '\n')
			t.write(','.join(self.pis.astype(str)) + '\n')
			t.write(str(np.sum(self.gamma_var*self.deltas*self.pis*self.KK)) + '\n')
			t.write(str(self.resid_var) + '\n')
			t.flush()

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
			for class_name in self.classes:
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


	def update_gamma_var(self, v0=0.0, s_sq=0.0, cc=1e-8):
		# First get middle gammas
		all_gammas_sq = []
		for window_name in self.window_names:
			window_sq_gammas = np.square(self.gamma[window_name])/self.deltas[self.class_membership[window_name]]
			#all_gammas_sq.append(window_sq_gammas[self.class_membership[window_name] != 0])
			all_gammas_sq.append(window_sq_gammas)
		all_gammas_sq = np.hstack(all_gammas_sq)

		vv = len(all_gammas_sq) + v0
		tau_sq = np.sum(all_gammas_sq) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		self.gamma_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_gamma(self):
		for window_name in self.window_names:
			gamma_vec, gwas_beta_resid_vec, class_membership_vec = update_gamma_from_single_window(np.load(self.window_info[window_name]['Q_file']), self.gamma[window_name], self.beta_resid[window_name], self.gamma_var, self.N_gwas, self.resid_var, self.pis, self.deltas)
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


		# Initialize variance parameters
		self.gamma_var = 1e-5

		# Initialize residual variance
		self.resid_var = 1.0

		# Number of categories
		self.C = len(self.deltas)

		# Self class names
		self.classes = np.arange(self.C)

		# Initialize pis
		self.pis = np.ones(self.C)/self.C


		# Initialize hyperparameter on pis
		self.alpha_0 = np.ones(self.C)


		# Keep track of sampled gamma_vars
		self.sampled_h2 = []
		self.sampled_resid_var = []


