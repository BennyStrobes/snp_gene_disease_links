import numpy as np
import os
import sys
import pdb
import gzip



def init_pis():
	pis = np.zeros(10)
	pis[0] = .314
	pis[1] = .1279
	pis[2] = .0879
	pis[3] = .0879
	pis[4] = .0879
	pis[5] = .0098
	pis[6] = .0098
	pis[7] = .0098
	pis[8] = .0098
	pis[9] = .0098
	pis = pis/np.sum(pis)
	return pis


def create_mapping_from_gene_name_to_gene_info(gene_annotation_file, chrom_num):
	f = open(gene_annotation_file)
	mapping = {}
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		if len(data) != 9:
			print('assumption eroror')
			pdb.set_trace()
		if data[2] != 'gene':
			continue
		ensamble_id = 'null'
		gene_type = 'null'
		gene_info = data[8].split(';')
		for info in gene_info:
			if info.startswith('gene_id'):
				ensamble_id = info.split('"')[1]
			elif info.startswith(' gene_type'):
				gene_type = info.split('"')[1]
		if ensamble_id == 'null' or gene_type == 'null':
			print('assumption eroror')
			pdb.set_trace()
		gene_chrom_num = data[0]
		if gene_chrom_num != 'chr' + chrom_num:
			continue
		gene_strand = data[6]
		if float(data[3]) > float(data[4]):
			print('assumption erroror')
			pdb.set_trace()
		if gene_strand == '+':
			tss = data[3]
		elif gene_strand == '-':
			tss = data[4]
		else:
			print('assumption error')


		# Add to info
		if ensamble_id not in mapping:
			mapping[ensamble_id] = (gene_type, gene_chrom_num, gene_strand, tss)
		else:
			if mapping[ensamble_id][0] != gene_type:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][1] != gene_chrom_num:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][2] != gene_strand:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][3] != tss:
				print('assumption eroror')
				pdb.set_trace()
	f.close()

	genes = []
	genes_tss = []
	for gene_name in [*mapping]:
		gene_info = mapping[gene_name]
		if gene_info[0] != 'protein_coding':
			continue
		if gene_info[1] != 'chr' + chrom_num:
			print('assumption eroror')
			pdb.set_trace()
		genes.append(gene_name)
		genes_tss.append(gene_info[3])

	genes = np.asarray(genes)
	gene_tss = np.asarray(genes_tss).astype(int)

	return genes, gene_tss


def closest_indexes(variant_pos, chrom_gene_positions, k=10):
	chrom_gene_positions = np.asarray(chrom_gene_positions)
	distances = np.abs(chrom_gene_positions - variant_pos)
	# Get indices of k smallest distances (unordered)
	closest_idx = np.argpartition(distances, k)[:k]
	# Sort those k by actual distance
	closest_idx = closest_idx[np.argsort(distances[closest_idx])]
	return closest_idx



def load_fine_mapped_variants(fm_summary_file):
	"""Read the PIP-filtered fine-mapping file. The bp column in this file is b37
	(confirmed 9/4/26: rs139913 bp=40713861 vs hg38 40317857), so positions here are
	IGNORED; hg38 positions come from the baselineLD annots via rsid instead."""
	variants = []
	f = open(fm_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if len(data) != 9:
			print('assumption eororr')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_chrom_num = data[1]
		rsid = data[3]
		pip = float(data[4])
		if pip < .01:
			print('assumpriororor')
			pdb.set_trace()
		variants.append((line_chrom_num, rsid, pip))
	f.close()
	return variants


def extract_hg38_variant_positions(baselineLD_anno_dir, needed_rsids):
	"""rsid -> (chrom, hg38 bp) from the same baselineLD annots used everywhere else in
	this pipeline. Only rsids in needed_rsids are kept so the mapping stays small."""
	mapping = {}
	for chrom_num in range(1,23):
		anno_file = baselineLD_anno_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'
		f = gzip.open(anno_file)
		head_count = 0
		for line in f:
			line = line.decode('utf-8').rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[2]
			if rsid not in needed_rsids:
				continue
			if rsid in mapping:
				continue
			mapping[rsid] = (data[0], int(data[1]))
		f.close()
	return mapping


######################
# command line args
######################
trait_name = sys.argv[1]
gene_annotation_file = sys.argv[2]
variant_fine_mapping_dir = sys.argv[3]
additive_sgl_output_file = sys.argv[4]
baselineLD_anno_dir = sys.argv[5]

# Pis
pis = init_pis()


# Extract gene info
gene_scores = {}
chrom_genes = {}
for chrom_num in range(1,23):
	gene_ensamble_ids, gene_positions = create_mapping_from_gene_name_to_gene_info(gene_annotation_file, str(chrom_num))
	for gene in gene_ensamble_ids:
		if gene in gene_scores:
			print('assumption oeroror')
			pdb.set_trace()
		gene_scores[gene] = 0.0
	chrom_genes[str(chrom_num)] = (gene_ensamble_ids, gene_positions)


# Load fine-mapped variants (the file's bp column is b37 and is ignored; see load_fine_mapped_variants)
fm_summary_file = variant_fine_mapping_dir + trait_name.split('460K.')[1] + '/' + 'PIP.postprocess.b37_filtered.txt'
fm_variants = load_fine_mapped_variants(fm_summary_file)

# hg38 position of each fine-mapped variant, keyed by rsid
rsid_to_hg38 = extract_hg38_variant_positions(baselineLD_anno_dir, set([var[1] for var in fm_variants]))

# Now loop through fine-mapped variants and update gene scores per each variant
n_missing = 0
n_chrom_mismatch = 0
for line_chrom_num, rsid, pip in fm_variants:
	if rsid not in rsid_to_hg38:
		n_missing = n_missing + 1
		continue
	anno_chrom, variant_pos = rsid_to_hg38[rsid]
	if anno_chrom != line_chrom_num:
		n_chrom_mismatch = n_chrom_mismatch + 1
		continue

	chrom_ensamble_ids, chrom_gene_positions = chrom_genes[str(line_chrom_num)]

	gene_indices = closest_indexes(variant_pos, chrom_gene_positions, k=10)

	k_closest_genes = chrom_ensamble_ids[gene_indices]

	line_gene_scores = pip*pis

	for gg, gene_name in enumerate(k_closest_genes):
		gene_scores[gene_name] = gene_scores[gene_name] + line_gene_scores[gg]

print('fine-mapped variants: ' + str(len(fm_variants)) + ' | no hg38 position in baselineLD (skipped): ' + str(n_missing) + ' | chromosome mismatch (skipped): ' + str(n_chrom_mismatch))



# Print to output file
t = open(additive_sgl_output_file,'w')
t.write('gene_name\tgene_score\n')
for gene_name in [*gene_scores]:
	t.write(gene_name + '\t' + str(gene_scores[gene_name]) + '\n')
t.close()
print(additive_sgl_output_file)


