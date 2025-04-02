import numpy as np
import os
import sys
import pdb








def tmp_loading(filer):
	gene_names_arr = []
	gene_name_to_integer = {}
	head_count = 0
	f = open(filer)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names_arr.append(data[0])
		gene_name_to_integer[data[0]] = int(data[1])
	return np.asarray(gene_names_arr), gene_name_to_integer


def convert_gene_names_mat_to_integer_mat(gene_names_mat, gene_name_to_integer):
	integer_mat = np.zeros(gene_names_mat.shape).astype(int)
	for ii in range(gene_names_mat.shape[0]):
		for jj in range(gene_names_mat.shape[1]):
			gene_name = gene_names_mat[ii,jj]
			integer_mat[ii,jj] = gene_name_to_integer[gene_name]

	return integer_mat

######################
# Command line args
######################
preorganized_snp_gene_annotation_dir = sys.argv[1]


input_file = preorganized_snp_gene_annotation_dir + 'window_LD_summary_with_snp_gene_anno.txt'

# Create mapping from gene name to integer
gene_name_to_integer = {}
gene_names_arr = []
counter = 0

head_count = 0
f = open(input_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene_names_file = data[7]
	gene_names_mat = np.load(gene_names_file)

	for ii in range(gene_names_mat.shape[0]):
		for jj in range(gene_names_mat.shape[1]):
			gene_name = gene_names_mat[ii,jj]
			if gene_name in gene_name_to_integer:
				continue
			gene_name_to_integer[gene_name] = counter
			gene_names_arr.append(gene_name)
			counter = counter + 1
f.close()
gene_names_arr = np.asarray(gene_names_arr)

# Print gene names arr to output
t = open(preorganized_snp_gene_annotation_dir + 'gene_name_to_integer_mapping.txt','w')
t.write('gene_name\tgene_integer\n')
for ii, gene_name in enumerate(gene_names_arr):
	t.write(gene_name + '\t' + str(ii) + '\n')
	if ii != gene_name_to_integer[gene_name]:
		print('assumption eroror')
		pdb.set_trace()
t.close()

#gene_names_arr, gene_name_to_integer = tmp_loading(preorganized_snp_gene_annotation_dir + 'gene_name_to_integer_mapping.txt')




input_file = preorganized_snp_gene_annotation_dir + 'window_LD_summary_with_snp_gene_anno.txt'
output_file = preorganized_snp_gene_annotation_dir + 'window_LD_summary_with_snp_gene_anno_v2.txt'


head_count = 0
f = open(input_file)
t = open(output_file,'w')


head_count = 0
f = open(input_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	print(data[0])
	gene_names_file = data[7]
	gene_names_mat = np.load(gene_names_file)

	gene_names_integer_mat = convert_gene_names_mat_to_integer_mat(gene_names_mat, gene_name_to_integer)
	gene_names_integer_file = gene_names_file.split('genes.npy')[0] + 'integer_genes.npy'
	np.save(gene_names_integer_file, gene_names_integer_mat)

	data = np.asarray(data)

	t.write('\t'.join(data[:7]) + '\t' + gene_names_integer_file + '\t' + data[8] + '\n')



f.close()
t.close()


