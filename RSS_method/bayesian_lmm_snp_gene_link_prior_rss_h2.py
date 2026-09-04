import numpy as np 
import pdb
import gzip
from scipy.stats import invgamma
import time
from numba import njit

@njit(cache=True)
def fast_log_sum_exp_vector(aa):
	a_max = np.max(aa)                     # scalar
	shifted = aa - a_max                   # avoid overflow
	sum_exp = np.sum(np.exp(shifted))    # scalar
	return a_max + np.log(sum_exp)

@njit(cache=True)
def fast_categorical_sample(probs):
	cum_probs = np.cumsum(probs)
	rr = np.random.rand()
	return np.searchsorted(cum_probs, rr)


def get_rsids_in_window(rsid_file):
	"""Read one window's snp ids. Called only when writing final output, so the ~5k strings
	are transient rather than 9.9M strings resident for the whole chain."""
	rsids = []
	f = open(rsid_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsids.append(line.rstrip().split('\t')[0])
	f.close()
	return rsids


def print_snp_gene_priors_output_file(pis, snp_gene_priors_output_file, has_null=False):
	"""Under the with-null method version pis has C+1 entries and entry 0 is the
	genome-wide background (null) class; it is written as prior_name 'null'."""
	t = open(snp_gene_priors_output_file,'w')
	t.write('prior_name\tprior_mean\tprior_mean_lb\tprior_mean_ub\n')
	for index, pi in enumerate(pis):
		if has_null:
			prior_name = 'null' if index == 0 else 'component_' + str(index - 1)
		else:
			prior_name = 'component_' + str(index)
		t.write(prior_name + '\t' + str(pi) + '\t' + str(pi) + '\t' + str(pi) + '\n')
	t.close()
	return


def print_gene_scores_to_output_file(sampled_gene_scores, ordered_genes, gene_score_output_file):
	"""Posterior summary of the per-gene score across the retained gibbs samples.

	The score is gamma_var[g]: the AVERAGE per-snp heritability among the snps linked to
	gene g. Its inverse-gamma conditional has mean ~ sum(gamma^2 over linked snps)/n_linked,
	so it does not grow with how many snps the gene captured. See the total-h2 output for
	the quantity that does.
	"""
	t = open(gene_score_output_file,'w')
	t.write('gene_name\tgene_score\tgene_score_lb\tgene_score_ub\tgene_score_variance\tgene_z_score\n')
	for g_index, gene_name in enumerate(ordered_genes):
		vec = sampled_gene_scores[:, g_index]
		# genes that never had a snp linked to them hold their initial value in every sample
		if len(np.unique(vec)) == 1:
			continue
		meaner = np.mean(vec)
		lower_bound, upper_bound = np.percentile(vec, [2.5, 97.5])
		variance = np.var(vec)
		t.write(gene_name + '\t' + str(meaner) + '\t' + str(lower_bound) + '\t' + str(upper_bound) + '\t' + str(variance) + '\t' + str(meaner/np.sqrt(variance)) + '\n')
	t.close()
	return


def print_gene_total_h2_to_output_file(sampled_gene_total_h2, mean_gene_n_snps, ordered_genes, gene_total_h2_output_file):
	"""Posterior summary of the TOTAL heritability linked to each gene.

	Where gene_score is the AVERAGE per-snp heritability among the snps linked to a gene,
	this is the SUM of gamma^2 over those snps -- the score scaled by how many snps the gene
	actually captured. A gene can score highly on one and not the other: a single strongly
	linked snp gives a high average but little total, and many weakly linked snps the reverse.
	"""
	t = open(gene_total_h2_output_file,'w')
	t.write('gene_name\tgene_total_h2\tgene_total_h2_lb\tgene_total_h2_ub\tgene_total_h2_variance\tgene_total_h2_z_score\tmean_n_linked_snps\n')
	for g_index, gene_name in enumerate(ordered_genes):
		vec = sampled_gene_total_h2[:, g_index]
		# genes that never captured a snp contribute no heritability in any sample
		if np.max(vec) == 0.0:
			continue
		meaner = np.mean(vec)
		lower_bound, upper_bound = np.percentile(vec, [2.5, 97.5])
		variance = np.var(vec)
		if variance > 0.0:
			zed = meaner/np.sqrt(variance)
		else:
			zed = np.nan
		t.write(gene_name + '\t' + str(meaner) + '\t' + str(lower_bound) + '\t' + str(upper_bound) + '\t' + str(variance) + '\t' + str(zed) + '\t' + str(mean_gene_n_snps[g_index]) + '\n')
	t.close()
	return


def print_snp_gene_link_output_file(window_names, snp_gene_class_counts, n_posterior_samples, window_info, ordered_genes, snp_gene_link_output_file):
	"""Posterior snp-gene assignment, written ONCE at the end of the run.

	Each snp's class assignment is tallied across the retained samples. Every candidate
	gene (all C ranks, closest first) is reported with its posterior probability -- one
	line per snp-gene pair, so the modal gene is the max-probability row per snp.
	Under the with-null method version counts has C+1 columns; column 0 is the background
	(no linked gene) class and is written with gene_rank -1 / gene_name null.
	Gzipped because this is C lines per snp genome-wide.
	"""
	t = gzip.open(snp_gene_link_output_file, 'wt')
	t.write('snp_id\tgene_rank\tgene_name\tgene_integer\tposterior_probability\n')
	for window_name in window_names:
		counts = snp_gene_class_counts[window_name]
		probs = counts/n_posterior_samples
		window_snp_gene_mat = np.load(window_info[window_name]['snp_gene_names_file'])
		window_rsids = get_rsids_in_window(window_info[window_name]['rsid_file'])
		CC = counts.shape[1]
		has_null = CC == window_snp_gene_mat.shape[1] + 1
		for snp_index, snp_name in enumerate(window_rsids):
			for snp_class in range(CC):
				if has_null and snp_class == 0:
					t.write(snp_name + '\t-1\tnull\t-1\t' + str(probs[snp_index, snp_class]) + '\n')
					continue
				gene_rank = snp_class - 1 if has_null else snp_class
				gene_index = window_snp_gene_mat[snp_index, gene_rank]
				t.write(snp_name + '\t' + str(gene_rank) + '\t' + ordered_genes[gene_index] + '\t' + str(gene_index) + '\t' + str(probs[snp_index, snp_class]) + '\n')
	t.close()
	return


@njit(cache=True)
def update_gamma_from_single_window(QQ_T, qq_sq, gamma_vec, gwas_beta_resid, gamma_var, N_gwas, resid_var, log_pis, snp_gene_link_mat):
	"""One gibbs sweep over the reference snps of a single window.

	QQ_T is Q TRANSPOSED: (n_reference_snps X n_whitened_components). Row j is snp j's
	loading vector and is contiguous, which is ~2.4x faster than striding down the columns
	of Q. qq_sq holds the squared column norms of Q; Q never changes across gibbs
	iterations so the caller computes it once.

	gamma_vec and gwas_beta_resid are updated IN PLACE and also returned.
	"""
	KK = len(gamma_vec)
	CC = snp_gene_link_mat.shape[1]
	NC = QQ_T.shape[1]

	class_membership_vec = np.zeros(KK, dtype=np.int32)

	posterior_var = np.zeros(CC)
	posterior_mean = np.zeros(CC)
	log_like = np.zeros(CC)

	scaler = N_gwas/resid_var

	for snp_index in np.random.permutation(KK):
		q_row = QQ_T[snp_index, :]
		old_gamma = gamma_vec[snp_index]

		# Re-include this snp's current effect, and accumulate resid.q in the same pass
		dotter = 0.0
		for cc in range(NC):
			gwas_beta_resid[cc] = gwas_beta_resid[cc] + q_row[cc]*old_gamma
			dotter = dotter + gwas_beta_resid[cc]*q_row[cc]

		posterior_var_term1 = qq_sq[snp_index]*scaler

		for cc in range(CC):
			gene_variance = gamma_var[snp_gene_link_mat[snp_index, cc]]
			pv = 1.0/(posterior_var_term1 + (1.0/gene_variance))
			pm = pv*scaler*dotter
			posterior_var[cc] = pv
			posterior_mean[cc] = pm
			if pv == 0.0:
				log_like[cc] = log_pis[cc]
			else:
				log_like[cc] = log_pis[cc] - (.5*np.log(gene_variance/pv)) + (.5*((pm*pm)/pv))

		probs = np.exp(log_like - fast_log_sum_exp_vector(log_like))

		class_p = fast_categorical_sample(probs)
		class_membership_vec[snp_index] = np.int32(class_p)

		new_gamma = np.random.normal(posterior_mean[class_p], np.sqrt(posterior_var[class_p]))
		gamma_vec[snp_index] = new_gamma

		# Remove the updated effect
		for cc in range(NC):
			gwas_beta_resid[cc] = gwas_beta_resid[cc] - q_row[cc]*new_gamma

	return gamma_vec, gwas_beta_resid, class_membership_vec



@njit(cache=True)
def update_gamma_from_single_window_with_null(QQ_T, qq_sq, gamma_vec, gwas_beta_resid, gamma_var, background_var, N_gwas, resid_var, log_pis, snp_gene_link_mat):
	"""One gibbs sweep of a single window under the with-null method version.

	Class 0 is a genome-wide background slab with variance background_var; classes 1..C
	are the C candidate genes (rank c-1 in snp_gene_link_mat). The class update is
	collapsed (gamma integrated out) exactly as in the original function -- the
	background is simply one more variance component in the categorical.
	"""
	KK = len(gamma_vec)
	CC = snp_gene_link_mat.shape[1]
	NC = QQ_T.shape[1]
	n_classes = CC + 1

	class_membership_vec = np.zeros(KK, dtype=np.int32)

	posterior_var = np.zeros(n_classes)
	posterior_mean = np.zeros(n_classes)
	log_like = np.zeros(n_classes)

	scaler = N_gwas/resid_var

	for snp_index in np.random.permutation(KK):
		q_row = QQ_T[snp_index, :]
		old_gamma = gamma_vec[snp_index]

		# Re-include this snp's current effect, and accumulate resid.q in the same pass
		dotter = 0.0
		for cc in range(NC):
			gwas_beta_resid[cc] = gwas_beta_resid[cc] + q_row[cc]*old_gamma
			dotter = dotter + gwas_beta_resid[cc]*q_row[cc]

		posterior_var_term1 = qq_sq[snp_index]*scaler

		for cc in range(n_classes):
			if cc == 0:
				gene_variance = background_var
			else:
				gene_variance = gamma_var[snp_gene_link_mat[snp_index, cc-1]]
			pv = 1.0/(posterior_var_term1 + (1.0/gene_variance))
			pm = pv*scaler*dotter
			posterior_var[cc] = pv
			posterior_mean[cc] = pm
			if pv == 0.0:
				log_like[cc] = log_pis[cc]
			else:
				log_like[cc] = log_pis[cc] - (.5*np.log(gene_variance/pv)) + (.5*((pm*pm)/pv))

		probs = np.exp(log_like - fast_log_sum_exp_vector(log_like))

		class_p = fast_categorical_sample(probs)
		class_membership_vec[snp_index] = np.int32(class_p)

		new_gamma = np.random.normal(posterior_mean[class_p], np.sqrt(posterior_var[class_p]))
		gamma_vec[snp_index] = new_gamma

		# Remove the updated effect
		for cc in range(NC):
			gwas_beta_resid[cc] = gwas_beta_resid[cc] - q_row[cc]*new_gamma

	return gamma_vec, gwas_beta_resid, class_membership_vec


class Bayesian_LMM_RSS_h2_inference(object):
	def __init__(self, window_names, window_info, n_gwas_individuals, ordered_genes, output_stem,inv_gamma_alpha=1e-10, inv_gamma_beta=1e-10, cross_gene_hyperparm=False,method_version='snp_gene_component_fixed_to_smart_init'):
		self.N_gwas = n_gwas_individuals
		self.window_names = window_names
		self.window_info = window_info
		self.ordered_genes = np.asarray(ordered_genes)
		self.output_stem = output_stem
		self.cross_gene_hyperparm = cross_gene_hyperparm

		self.inv_gamma_alpha = inv_gamma_alpha
		self.inv_gamma_beta = inv_gamma_beta
		supported_method_versions = ['snp_gene_component_fixed_to_smart_init', 'snp_gene_component_with_null']
		if method_version not in supported_method_versions:
			print('assumption eroror: unsupported method_version ' + str(method_version) + '; supported: ' + ','.join(supported_method_versions))
			pdb.set_trace()
		if method_version == 'snp_gene_component_with_null' and not cross_gene_hyperparm:
			# the sparse model redraws empty genes' slab variances from their prior each sweep,
			# which is only proper/stable under the learned cross-gene hyperprior
			print('assumption eroror: snp_gene_component_with_null requires prior_choice=inverse_gamma_cross_gene_prior_*')
			pdb.set_trace()
		self.method_version = method_version

		# Number of snps
		self.KK = 0.0
		for window in window_names:
			self.KK = self.KK + self.window_info[window]['n_snps']
		
		return



	def fit(self, total_iterations=10000, burn_in_iterations=5000, update_resid_var_bool=True, thin_iterations=5):
		""" Fit the model.
		"""
		# Initialize model params
		self.initialize_variables()

		# Preallocate the posterior sample arrays. Appending to a list and calling np.asarray
		# at the end would transiently hold both copies; float32 is ample for reporting means,
		# percentiles and variances.
		n_expected_samples = 0
		for itera in range(total_iterations):
			if itera > burn_in_iterations and np.mod(itera + 1, thin_iterations) == 0:
				n_expected_samples = n_expected_samples + 1
		self.sampled_gene_scores = np.zeros((n_expected_samples, self.G), dtype=np.float32)
		self.sampled_gene_total_h2 = np.zeros((n_expected_samples, self.G), dtype=np.float32)
		self.sampled_resid_var = np.zeros(n_expected_samples)

		# Keep track of iterations
		self.itera = 0

		#print(self.output_log_file)
		t = open(self.output_stem + 'log.txt','w')
		print(self.output_stem + 'log.txt')

		# Iterative Gibbs sampling algorithm
		t1 = time.time()
		for itera in range(total_iterations):
			print(itera)
			
			# Update gamma
			self.update_gamma()
			t2 = time.time()
			print(t2-t1)
			t1 = time.time()

			# pis are held FIXED at their smart-init values under the original method version.
			# Under snp_gene_component_with_null only the null probability pi0 and the
			# background variance are learned genome-wide each sweep; conditional on non-null
			# the split over the 10 gene ranks stays fixed at the smart-init values.
	
			if self.method_version == 'snp_gene_component_with_null':
				self.update_pi()
				self.update_background_var()

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
			if self.method_version == 'snp_gene_component_with_null':
				t.write(str(self.background_var) + '\n')
			t.flush()
			print(','.join(self.pis.astype(str)))


			# Retain a thinned posterior sample instead of dumping a full set of output files
			# every few iterations. Summaries are written once, after the chain finishes.
			if itera > burn_in_iterations and np.mod(self.itera, thin_iterations) == 0:
				self.record_posterior_sample()

			'''
			if np.mod(self.itera,100) == 0.0:
				print(self.gamma_var*self.KK)
			'''
			#t.write('Iteration: ' + str(itera) + '\n')
			#t.write(','.join(self.pis.astype(str)) + '\n')
			#t.write(str(np.sum(self.gamma_var*self.deltas*self.pis*self.KK)) + '\n')
			#t.write(str(self.resid_var) + '\n')
			#t.flush()


		t.close()

		###########################
		# Posterior summaries, written ONCE
		###########################
		# a short-fall here would leave preallocated zeros polluting every posterior mean
		if self.n_posterior_samples != self.sampled_gene_scores.shape[0]:
			print('assumption eroror: retained ' + str(self.n_posterior_samples) + ' posterior samples but preallocated ' + str(self.sampled_gene_scores.shape[0]))
			pdb.set_trace()
		if self.n_posterior_samples < 2:
			print('assumption eroror: only ' + str(self.n_posterior_samples) + ' posterior samples retained; check total_iterations / burn_in_iterations / thin_iterations')
			pdb.set_trace()

		print('retained ' + str(self.n_posterior_samples) + ' posterior samples (burn in ' + str(burn_in_iterations) + ', thin ' + str(thin_iterations) + ')')

		print_gene_scores_to_output_file(self.sampled_gene_scores, self.ordered_genes, self.output_stem + 'gene_score_averaged.txt')
		print_gene_total_h2_to_output_file(self.sampled_gene_total_h2, self.sum_gene_n_snps/self.n_posterior_samples, self.ordered_genes, self.output_stem + 'gene_total_h2_averaged.txt')
		print_snp_gene_priors_output_file(self.pis, self.output_stem + 'snp_gene_priors_averaged.txt', has_null=(self.method_version == 'snp_gene_component_with_null'))
		print_snp_gene_link_output_file(self.window_names, self.snp_gene_class_counts, self.n_posterior_samples, self.window_info, self.ordered_genes, self.output_stem + 'snp_gene_links.txt.gz')

		return

	def record_posterior_sample(self):
		ii = self.n_posterior_samples
		self.sampled_gene_scores[ii, :] = self.gamma_var
		self.sampled_gene_total_h2[ii, :] = self.gene_total_h2
		self.sampled_resid_var[ii] = self.resid_var

		# only the MEAN linked-snp count is reported, so a running sum suffices
		self.sum_gene_n_snps = self.sum_gene_n_snps + self.gene_n_snps

		# Tally each snp's gene assignment so the final links file can report a posterior
		# mode rather than one arbitrary draw
		for window_name in self.window_names:
			counts = self.snp_gene_class_counts[window_name]
			counts[np.arange(counts.shape[0]), self.class_membership[window_name]] += 1

		self.n_posterior_samples = self.n_posterior_samples + 1

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

	def update_pi(self):
		# Genome-wide conjugate beta-bernoulli update of the SINGLE null probability pi0.
		# Conditional on non-null, the split over the C gene ranks stays FIXED at the
		# smart-init values (matching the original method version and the additive baseline).
		n_null = 0.0
		n_gene = 0.0
		for window_name in self.window_names:
			cm = self.class_membership[window_name]
			n_null = n_null + np.sum(cm == 0)
			n_gene = n_gene + np.sum(cm > 0)
		pi0 = np.random.beta(self.pi0_beta_prior_a + n_null, self.pi0_beta_prior_b + n_gene)
		self.pis = np.hstack(([pi0], (1.0 - pi0)*self.smart_init_pis))
		self.log_pis = np.log(self.pis)
		return

	def update_background_var(self, cc=1e-10):
		# Conjugate inverse-gamma update of the genome-wide background (null-class) variance
		n_background = 0.0
		sum_sq = 0.0
		for window_name in self.window_names:
			background_snps = self.class_membership[window_name] == 0
			n_background = n_background + np.sum(background_snps)
			sum_sq = sum_sq + np.sum(np.square(self.gamma[window_name][background_snps]))
		if n_background < 1:
			# no snps currently in the background class; keep the current value
			return
		self.background_var = 1.0/np.random.gamma(shape=(n_background/2) + cc, scale=1.0/((sum_sq/2) + cc))
		return

	def update_gamma_var(self, v0=0.0, s_sq=0.0):
		# First get middle gammas
		sum_gamma_sq_vec = np.zeros(self.G)
		num_snps = np.zeros(self.G)
		with_null = self.method_version == 'snp_gene_component_with_null'
		for window_name in self.window_names:
			window_sq_gammas = np.square(self.gamma[window_name])
			window_class_membership = self.class_membership[window_name]
			window_snp_gene_mat = np.load(self.window_info[window_name]['snp_gene_names_file'])

			if with_null:
				# class 0 is the background component: those snps belong to no gene
				gene_snp_indices = np.where(window_class_membership > 0)[0]
				linked_genes = window_snp_gene_mat[gene_snp_indices, window_class_membership[gene_snp_indices] - 1]
				np.add.at(sum_gamma_sq_vec, linked_genes, window_sq_gammas[gene_snp_indices])
				np.add.at(num_snps, linked_genes, np.ones(len(gene_snp_indices)))
			else:
				linked_genes = window_snp_gene_mat[np.arange(window_snp_gene_mat.shape[0]),window_class_membership]
				np.add.at(sum_gamma_sq_vec, linked_genes, window_sq_gammas)
				np.add.at(num_snps, linked_genes, np.ones(len(window_sq_gammas)))




		vv = num_snps + v0
		tau_sq = sum_gamma_sq_vec + s_sq

		self.valid_gene_indices = num_snps >= 1


		param1 = (vv/2) + self.inv_gamma_alpha # inv_gamma_alpha == gene_gamma_shape
		param2 = 1.0/((tau_sq/2) + self.inv_gamma_beta) # inv_gamma_beta == self.gene_gamma_scale)


		self.gamma_var[self.valid_gene_indices] = 1.0/np.random.gamma(shape=param1[self.valid_gene_indices], scale=param2[self.valid_gene_indices])

		if with_null:
			# Under the sparse model many genes hold zero snps in a given sweep. Their
			# conditional posterior IS the (hierarchically shrunk) prior, so redraw from it
			# rather than freezing a stale value. Prior draws are proper here because the
			# with-null version requires the learned cross-gene hyperprior on inv_gamma_beta.
			n_empty = int(np.sum(~self.valid_gene_indices))
			if n_empty > 0:
				self.gamma_var[~self.valid_gene_indices] = 1.0/np.random.gamma(shape=self.inv_gamma_alpha, scale=1.0/self.inv_gamma_beta, size=n_empty)

		# Retained for the per-gene total-h2 output: the realized sum of squared effects
		# currently linked to each gene, and how many snps that was
		self.gene_total_h2 = sum_gamma_sq_vec
		self.gene_n_snps = num_snps

		return

	def update_gamma(self):
		for window_name in self.window_names:
			# Q is (n_components X n_reference_snps). The gibbs loop walks reference snps, so
			# transposing makes each snp's loading vector contiguous. Measured ~2.4x faster than
			# striding down columns, and the transpose costs less than the difference.
			QQ_T = np.ascontiguousarray(np.transpose(np.load(self.window_info[window_name]['Q_file'])), dtype=np.float64)

			# Q is constant across gibbs iterations, so its squared column norms are computed once
			if 'qq_sq' not in self.window_info[window_name]:
				self.window_info[window_name]['qq_sq'] = np.sum(np.square(QQ_T), axis=1)

			if self.method_version == 'snp_gene_component_with_null':
				gamma_vec, gwas_beta_resid_vec, class_membership_vec = update_gamma_from_single_window_with_null(QQ_T, self.window_info[window_name]['qq_sq'], self.gamma[window_name], self.beta_resid[window_name], self.gamma_var, self.background_var, self.N_gwas, self.resid_var, self.log_pis, np.load(self.window_info[window_name]['snp_gene_names_file']))
			else:
				gamma_vec, gwas_beta_resid_vec, class_membership_vec = update_gamma_from_single_window(QQ_T, self.window_info[window_name]['qq_sq'], self.gamma[window_name], self.beta_resid[window_name], self.gamma_var, self.N_gwas, self.resid_var, self.log_pis, np.load(self.window_info[window_name]['snp_gene_names_file']))
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
			self.class_membership[window_name] = np.zeros(self.window_info[window_name]['n_snps'], dtype=np.int32)


		random_snp_gene_names_mat = np.load(self.window_info[window_name]['snp_gene_names_file'])

		# One category per linked gene -- there is no null component in this method version
		self.C = random_snp_gene_names_mat.shape[1]

		# Number of genes
		self.G = len(self.ordered_genes)

		# Initialize variance parameters
		self.gamma_var = np.ones(self.G)*.1/1000000

		# Initialize residual variance
		self.resid_var = 1.0

		# Initialize pis
		# The smart init hardcodes a prior over the 10 nearest genes
		if self.C != 10:
			print('assumption eroror: smart init expects 10 linked genes per snp, got ' + str(self.C))
			pdb.set_trace()
		smart_init_pis = np.ones(self.C)/self.C
		if True:
			smart_init_pis[0] = .314
			smart_init_pis[1] = .1279
			smart_init_pis[2] = .0879
			smart_init_pis[3] = .0879
			smart_init_pis[4] = .0879
			smart_init_pis[5] = .0098
			smart_init_pis[6] = .0098
			smart_init_pis[7] = .0098
			smart_init_pis[8] = .0098
			smart_init_pis[9] = .0098
			smart_init_pis = smart_init_pis/np.sum(smart_init_pis)

		if self.method_version == 'snp_gene_component_with_null':
			# Class 0 is a genome-wide background component (diffuse per-snp heritability, no
			# linked gene); classes 1..C are the C candidate genes.
			self.n_classes = self.C + 1
			# Only pi0 (null vs any-gene) is learned; conditional on non-null the split over
			# the C gene ranks stays fixed at the smart-init values for the whole run.
			self.smart_init_pis = smart_init_pis
			pi0_init = 0.9
			self.pis = np.hstack(([pi0_init], (1.0 - pi0_init)*smart_init_pis))
			# weak beta prior centered on the init; genome-wide counts dominate it immediately
			self.pi0_beta_prior_a = 9.0
			self.pi0_beta_prior_b = 1.0
			# Initialization ORDER is what makes the background/slab roles identifiable: the
			# background starts at the average per-snp effect scale (mean beta^2 minus the
			# sampling-noise floor) and the gene slabs start 100x above it, so large effects
			# are claimed by genes from the first sweep. The conjugate variance updates move
			# down quickly but crawl up only slowly, so the background can never climb into
			# the slab role (the label-swapped mode is unreachable) while an overestimated
			# init self-corrects downward within a few hundred sweeps.
			sum_sq_beta = 0.0
			n_beta = 0.0
			for window_name in self.window_names:
				sum_sq_beta = sum_sq_beta + np.sum(np.square(self.window_info[window_name]['beta_pc']))
				n_beta = n_beta + len(self.window_info[window_name]['beta_pc'])
			self.background_var = max(sum_sq_beta/n_beta - 1.0/self.N_gwas, 1e-3/self.N_gwas)
			self.gamma_var = np.ones(self.G)*self.background_var*100.0
		else:
			self.n_classes = self.C
			# pis are fixed for the whole run under this method version
			self.pis = smart_init_pis

		self.log_pis = np.log(self.pis)


		# Posterior accumulators. These replace writing a full set of output files every few
		# gibbs iterations; summaries are computed from them once the chain finishes.
		self.n_posterior_samples = 0
		self.sum_gene_n_snps = np.zeros(self.G)

		# uint16 is exact up to 65535 retained samples and halves this array, which is the
		# largest thing the sampler holds
		self.snp_gene_class_counts = {}
		for window_name in self.window_names:
			self.snp_gene_class_counts[window_name] = np.zeros((self.window_info[window_name]['n_snps'], self.n_classes), dtype=np.uint16)

		# filled in by update_gamma_var
		self.gene_total_h2 = np.zeros(self.G)
		self.gene_n_snps = np.zeros(self.G)

		return


