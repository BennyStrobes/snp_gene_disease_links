# Synthetic-data validation of method_version='snp_gene_component_with_null'
# (and regression check of the original fixed-pis method). Run: python3 test_with_null_method_synthetic.py
#
# Simulates 4 independent windows (identity LD), 40 genes, with true pi0=0.85,
# background var 5e-6 (near the noise floor on purpose) and slab var 2e-4.
# Notes on expectations:
#  - pi0 is identified only up to an indistinguishability band (~0.75-0.96): sub-noise
#    slab draws are genuinely null-like, and near-noise background effects leak out.
#  - background_var initializes at the average per-snp effect scale and should settle
#    near the true background variance, well below the gene slabs (init ordering is
#    what makes the background/slab roles identifiable).
#  - under identity LD only the SUM of resid_var and N*pi0*background_var is well
#    identified for diffuse snps; the test asserts on that combined quantity.
import sys, os, gzip, tempfile, contextlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numba import njit
import bayesian_lmm_snp_gene_link_prior_rss_h2 as blmm

@njit(cache=False)
def _seed_numba(s):
	np.random.seed(s)

np.random.seed(1)
_seed_numba(7)

N_GWAS = 100000
N_WIN, N_SNP, N_GENE = 4, 300, 40
PI0_TRUE, BG_VAR_TRUE, SLAB_VAR_TRUE = 0.85, 5e-6, 2e-4

d = tempfile.mkdtemp()
genes = ['gene%d' % i for i in range(N_GENE)]
window_names, window_info, truth_null = [], {}, {}
for w in range(N_WIN):
	wn = 'w%d' % w
	window_names.append(wn)
	rsid_f = os.path.join(d, wn + '_rsids.txt')
	with open(rsid_f, 'w') as t:
		t.write('rsid\ta0\ta1\n')
		for s in range(N_SNP):
			t.write('rs_%d_%d\tA\tG\n' % (w, s))
	np.save(os.path.join(d, wn + '_sg.npy'), np.tile(np.arange(w*10, w*10+10), (N_SNP, 1)))
	np.save(os.path.join(d, wn + '_Q.npy'), np.eye(N_SNP))
	is_null = np.random.rand(N_SNP) < PI0_TRUE
	gamma_true = np.where(is_null, np.random.normal(0, np.sqrt(BG_VAR_TRUE), N_SNP),
	                      np.random.normal(0, np.sqrt(SLAB_VAR_TRUE), N_SNP))
	beta_pc = gamma_true + np.random.normal(0, np.sqrt(1.0/N_GWAS), N_SNP)
	truth_null[wn] = is_null
	window_info[wn] = {'n_snps': N_SNP, 'rsid_file': rsid_f, 'Q_file': os.path.join(d, wn + '_Q.npy'),
	                   'beta_pc': beta_pc, 'snp_gene_names_file': os.path.join(d, wn + '_sg.npy')}

# ---------------- with-null path ----------------
stem = os.path.join(d, 'new_')
mod = blmm.Bayesian_LMM_RSS_h2_inference(np.asarray(window_names), window_info, N_GWAS, genes, stem,
        inv_gamma_alpha=1.0, inv_gamma_beta=1e-7, cross_gene_hyperparm=True,
        method_version='snp_gene_component_with_null')
with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
	mod.fit(total_iterations=400, burn_in_iterations=200, update_resid_var_bool=True, thin_iterations=2)

pi0s, bgs = [], []
lines = [l for l in open(stem + 'log.txt') if not l.startswith('Iteration')]
for i in range(0, len(lines), 4):
	pi0s.append(float(lines[i].split(',')[0])); bgs.append(float(lines[i+3]))
pi0_mean, bg_mean = np.mean(pi0s[-200:]), np.mean(bgs[-200:])
print('posterior-mean pi0: %.3f (truth %.2f; identified band ~0.75-0.96)' % (pi0_mean, PI0_TRUE))
assert 0.72 < pi0_mean < 0.96
print('posterior-mean background_var: %.2e (truth %.0e; must stay << slab)' % (bg_mean, BG_VAR_TRUE))
assert 5e-7 < bg_mean < 5e-5
null_cm = np.hstack([mod.class_membership[wn][truth_null[wn]] for wn in window_names])
gene_cm = np.hstack([mod.class_membership[wn][~truth_null[wn]] for wn in window_names])
f_null, f_gene = np.mean(null_cm == 0), np.mean(gene_cm == 0)
print('frac in class0 | truly null: %.3f | truly linked: %.3f' % (f_null, f_gene))
assert 0.75 < f_null < 0.98 and (f_null - f_gene) > 0.25
# Under identity LD, background_var and resid_var/N sit on a near-ridge (only their sum
# is well identified for diffuse snps), so assert on the COMBINED diffuse variance:
# truth = 1.0 noise + N*pi0*bg = 1.43. Real LD separates the two better than this synthetic.
combined = mod.resid_var + N_GWAS*pi0_mean*bg_mean
print('resid_var: %.3f | combined diffuse variance: %.3f (truth ~1.43)' % (mod.resid_var, combined))
assert 0.3 < mod.resid_var < 2.2
assert 1.0 < combined < 1.9
rank0 = [0, 10, 20, 30]
print('mean slab var, rank-0 genes: %.2e (truth %.0e)' % (np.mean(mod.gamma_var[rank0]), SLAB_VAR_TRUE))
assert 3e-5 < np.mean(mod.gamma_var[rank0]) < 2e-3

th2 = {}
for i, line in enumerate(open(stem + 'gene_total_h2_averaged.txt')):
	if i == 0:
		continue
	p = line.split('\t'); th2[p[0]] = float(p[1])
top8 = sorted(th2, key=th2.get, reverse=True)[:8]
assert all(('gene%d' % g) in top8 for g in rank0), top8
print('all 4 rank-0 genes in top-8 by total h2  OK')

rows = open(stem + 'snp_gene_priors_averaged.txt').readlines()
assert rows[1].split('\t')[0] == 'null' and len(rows) == 12
print('priors file: null row + 10 component rows  OK')

null_prob, nrows = {}, 0
for i, line in enumerate(gzip.open(stem + 'snp_gene_links.txt.gz', 'rt')):
	if i == 0:
		continue
	p = line.rstrip().split('\t'); nrows += 1
	if p[1] == '-1':
		assert p[2] == 'null' and p[3] == '-1'
		null_prob[p[0]] = float(p[4])
assert nrows == N_WIN*N_SNP*11
np_null = np.mean([null_prob['rs_%d_%d' % (w, s)] for w in range(N_WIN) for s in range(N_SNP) if truth_null['w%d' % w][s]])
np_gene = np.mean([null_prob['rs_%d_%d' % (w, s)] for w in range(N_WIN) for s in range(N_SNP) if not truth_null['w%d' % w][s]])
print('links: 11 rows/snp; mean null prob | truly null %.3f, truly linked %.3f' % (np_null, np_gene))
assert (np_null - np_gene) > 0.25

# ---------------- old-path regression ----------------
smart = np.array([.314, .1279, .0879, .0879, .0879, .0098, .0098, .0098, .0098, .0098]); smart /= smart.sum()
stem_old = os.path.join(d, 'old_')
mod2 = blmm.Bayesian_LMM_RSS_h2_inference(np.asarray(window_names), window_info, N_GWAS, genes, stem_old,
        inv_gamma_alpha=1e-16, inv_gamma_beta=1e-16, cross_gene_hyperparm=False,
        method_version='snp_gene_component_fixed_to_smart_init')
with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
	mod2.fit(total_iterations=40, burn_in_iterations=20, update_resid_var_bool=True, thin_iterations=5)
assert np.allclose(mod2.pis, smart)
lines = gzip.open(stem_old + 'snp_gene_links.txt.gz', 'rt').readlines()
assert len(lines) == 1 + N_WIN*N_SNP*10
assert not any('\tnull\t' in l for l in lines[1:])
assert open(stem_old + 'snp_gene_priors_averaged.txt').readlines()[1].startswith('component_0')
print('old-path regression: pis fixed, 10 rows/snp, no null rows  OK')
print('ALL VALIDATION CHECKS PASSED')
