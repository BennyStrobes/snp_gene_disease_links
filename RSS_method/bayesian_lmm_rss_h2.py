import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time

def update_gamma_from_single_window(QQ, gamma_vec, gwas_beta_resid, gamma_var, N_gwas, resid_var, univariate=True):
	if univariate == False:
		KK = len(gamma_vec)
		SS_inv = ((LD*N_gwas) + np.eye(KK)/gamma_var)
		SS= np.linalg.inv(SS_inv)
		pdb.set_trace()
		# NOTE DONE NEET TO ACCOUNT QT
		posterior_mean = N_gwas*np.dot(SS,gwas_beta_resid)
		gamma_vec = np.random.multivariate_normal(mean=posterior_mean, cov=SS)

	else:
		# N-snps in window
		KK = len(gamma_vec)

		# Precompute posterior variance (same for all snps)
		qq_sq = np.sum(np.square(QQ),axis=0)
		posterior_var = 1.0/((qq_sq*N_gwas/resid_var) + (1.0/gamma_var))

		# Update each snp in turn
		for snp_index in np.random.permutation(range(KK)):
			# Re include current effect
			gwas_beta_resid = gwas_beta_resid + (QQ[:, snp_index]*gamma_vec[snp_index])

			# Compute posterior mean for this snp
			posterior_mean = posterior_var[snp_index]*(N_gwas/resid_var)*np.dot(gwas_beta_resid,QQ[:,snp_index])

			# Sample from posterior distribution
			gamma_vec[snp_index] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var[snp_index]))

			# Remove updated effect
			gwas_beta_resid = gwas_beta_resid - (QQ[:, snp_index]*gamma_vec[snp_index])

	return gamma_vec, gwas_beta_resid





class Bayesian_LMM_RSS_h2_inference(object):
	def __init__(self, window_names, window_info, n_gwas_individuals, output_log_file):
		self.N_gwas = n_gwas_individuals
		self.window_names = window_names
		self.window_info = window_info
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
			t.write(str(self.gamma_var*self.KK) + '\n')
			t.write(str(self.resid_var) + '\n')
			t.flush()

			if itera > burn_in_iterations:
				self.sampled_h2.append(self.gamma_var*self.KK)
				self.sampled_resid_var.append(self.resid_var)


		self.sampled_h2 = np.asarray(self.sampled_h2)
		self.sampled_resid_var = np.asarray(self.sampled_resid_var)

		t.close()
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
		all_gammas = []
		for window_name in self.window_names:
			window_gammas = self.gamma[window_name]
			all_gammas.append(window_gammas)
		all_gammas = np.hstack(all_gammas)

		vv = len(all_gammas) + v0
		tau_sq = np.sum(np.square(all_gammas)) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		self.gamma_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_gamma(self):
		for window_name in self.window_names:

			gamma_vec, gwas_beta_resid_vec = update_gamma_from_single_window(np.load(self.window_info[window_name]['Q_file']), self.gamma[window_name], self.beta_resid[window_name], self.gamma_var, self.N_gwas, self.resid_var)
			self.gamma[window_name] = gamma_vec
			self.beta_resid[window_name] = gwas_beta_resid_vec

		return


	def initialize_variables(self):
		# Initialize causal effcts
		self.gamma = {}
		self.beta_resid = {}

		for window_name in self.window_names:
			self.gamma[window_name] = np.zeros(self.window_info[window_name]['n_snps'])
			self.beta_resid[window_name] = np.copy(self.window_info[window_name]['beta_pc'])

		# Initialize variance parameters
		self.gamma_var = 1e-5

		# Initialize residual variance
		self.resid_var = 1.0

		# Keep track of sampled gamma_vars
		self.sampled_h2 = []
		self.sampled_resid_var = []