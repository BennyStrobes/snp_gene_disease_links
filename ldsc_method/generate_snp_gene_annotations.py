import numpy as np
import sys
import pdb
import os


def get_distance_bins(snp_gene_annotation_dir,n_distance_bins):
	distances = []
	for chrom_num in range(20,23):
		print(chrom_num)
		filer = snp_gene_annotation_dir + 'snp_gene_pairs_chr' + str(chrom_num) + '.txt'
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			distances.append(np.asarray(data[4].split(';')).astype(float))
		f.close()
	distances = np.sort(np.asarray(np.hstack(distances)))
	bins = np.array_split(distances,n_distance_bins)

	bin_ranges = []
	for bin_values in bins:
		bin_ranges.append((np.min(bin_values), np.max(bin_values)))
	bin_ranges = np.asarray(bin_ranges)

	return bin_ranges








#####################
# Command line args
######################
snp_gene_annotation_dir = sys.argv[1]


# Get distance bins
n_distance_bins = 10
distance_bins = get_distance_bins(snp_gene_annotation_dir,n_distance_bins)

# Open output annotation file
output_annotation_file = snp_gene_annotation_dir + 'snp_gene_pair_annotations.txt'
t = open(output_annotation_file,'w')
for chrom_num in range(1,23):
	print(chrom_num)
	f = open(snp_gene_annotation_dir + 'snp_gene_pairs_chr' + str(chrom_num) + '.txt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		distances = np.asarray(data[4].split(';')).astype(float)
		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3])

		for gene_index,gene_distance in enumerate(distances):
			gene_distance_feature = ((gene_distance >= distance_bins[:,0]) & (gene_distance <= distance_bins[:,1]))*1.0
			gene_rank_feature = np.zeros(len(distances))
			gene_rank_feature[gene_index] = 1.0
			gene_feature_vector = np.hstack((gene_distance_feature, gene_rank_feature))
			t.write('\t' + ';'.join(gene_feature_vector.astype(str)))
		t.write('\n')


	f.close()

t.close()