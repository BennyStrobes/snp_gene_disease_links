import numpy as np
import os
import sys
import pdb










###################
# Command line args
###################
preorganized_snp_gene_annotation_dir = sys.argv[1]


output_file = preorganized_snp_gene_annotation_dir + 'window_LD_summary_with_snp_gene_anno_v2.txt'


head_count = 0


head_count = 0
gene_counts = np.zeros(17970)
max_gene = 0
f = open(output_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene_names_file = data[7]
	gene_names_mat = np.load(gene_names_file)

	closest_genes = gene_names_mat[:,0]

	for closest_gene_to_snp in closest_genes:
		gene_counts[closest_gene_to_snp] = gene_counts[closest_gene_to_snp] + 1

	max_gene = np.max([np.max(gene_names_mat), max_gene])

f.close()

