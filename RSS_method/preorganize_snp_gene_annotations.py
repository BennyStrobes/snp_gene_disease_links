import numpy as np
import os
import sys
import pdb







def create_mapping_from_snps_to_genes_and_anno_matrix(filer, chrom_num):
	# The reference (input) snp set is no longer restricted to hapmap3, so the annotation
	# file can legitimately repeat an rsid (eg. multi-allelic sites sharing an rsid in the
	# baselineLD annot file). A pdb trap here would hang a batch job, so count instead.
	dicti = {}
	n_identical_dups = 0
	n_conflicting_dups = 0
	f = open(filer)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip variants not on chromosome
		if data[1] != chrom_num:
			continue

		# Extract_relevent fields
		variant_id = data[0]
		gene_names = np.asarray(data[3].split(';'))
		snp_gene_feature_matrix = []
		for ele in data[4:]:
			snp_gene_feature_matrix.append(np.asarray(ele.split(';')).astype(float))
		snp_gene_feature_matrix = np.asarray(snp_gene_feature_matrix)

		if variant_id in dicti:
			# Keep the first occurrence either way; only shout when the repeat disagrees.
			if np.array_equal(dicti[variant_id][0], gene_names):
				n_identical_dups = n_identical_dups + 1
			else:
				n_conflicting_dups = n_conflicting_dups + 1
				if n_conflicting_dups <= 10:
					print('WARNING: conflicting duplicate annotation for ' + variant_id + ' (keeping first)')
			continue

		dicti[variant_id] = (gene_names, snp_gene_feature_matrix)

	f.close()

	if n_identical_dups > 0:
		print('chr' + chrom_num + ': collapsed ' + str(n_identical_dups) + ' identical duplicate rsids')
	if n_conflicting_dups > 0:
		print('##################################################################')
		print('chr' + chrom_num + ': ' + str(n_conflicting_dups) + ' CONFLICTING duplicate rsids (kept first occurrence)')
		print('##################################################################')

	return dicti





def extract_window_rsids(rsid_file):
	f = open(rsid_file)
	rsids = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsids.append(data[0])
	f.close()
	return np.asarray(rsids)


#######################
# Command line args
#######################
snp_gene_annotation_dir = sys.argv[1]
processed_ld_data_dir = sys.argv[2]
preorganized_snp_gene_annotation_dir = sys.argv[3]

# Must match the snp_set used to build the LD matrices in process_LD_data.sh
snp_set='assymetric_snps' # alt: hampmap3_snps


output_file = preorganized_snp_gene_annotation_dir + 'window_LD_summary_with_snp_gene_anno.txt'
t = open(output_file,'w')

# Loop through chromosomes
for chrom_num in range(1,23):
	# Create mapping from snp to genes and anno matrix
	snp_to_genes_and_anno_matrix = create_mapping_from_snps_to_genes_and_anno_matrix(snp_gene_annotation_dir + 'snp_gene_pair_annotations.txt', str(chrom_num))
	

	ld_summary_file = processed_ld_data_dir + 'window_LD_summary_' + snp_set + '_chr' + str(chrom_num) + '.txt'

	f = open(ld_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			if chrom_num == 1:
				t.write(line + '\t' + 'linked_genes_file\tsnp_gene_anno_file\n')
			continue
		# Extract relevent fields
		window_name = data[0]
		# Column 4 is the REFERENCE (input) snp file. Those snps index the COLUMNS of Q_mat,
		# ie. they are the snps whose causal effects (gamma) are estimated and linked to genes,
		# so they are what the rows of the saved gene annotations must line up with.
		# Column 5 is the regression (output) snp file -- the snps carrying gwas z-scores, which
		# index the ROWS of Q_mat. Under asymmetric LD these two sets differ.
		rsid_file = data[4]
		print(window_name)

		# Extract window rsids
		window_rsids = extract_window_rsids(rsid_file)

		# Extract data for window
		nearby_genes = []
		variant_gene_anno = []

		for rsid in window_rsids:
			if rsid not in snp_to_genes_and_anno_matrix:
				print('assumption eroror: reference snp ' + rsid + ' in window ' + window_name + ' has no snp-gene annotation.')
				print('The snp-gene annotation file must cover the full reference snp set, not just hapmap3 snps.')
				pdb.set_trace()
			gene_info = snp_to_genes_and_anno_matrix[rsid]
			nearby_genes.append(gene_info[0])
			variant_gene_anno.append(gene_info[1])

		if len(nearby_genes) == 0:
			print('assumption eroror: window ' + window_name + ' has no reference snps')
			pdb.set_trace()

		nearby_genes = np.asarray(nearby_genes) # N_snpsX10
		variant_gene_anno = np.asarray(variant_gene_anno)

		# Save to output
		linked_genes_file = preorganized_snp_gene_annotation_dir + window_name + '_linked_genes.npy'
		np.save(linked_genes_file, nearby_genes)
		snp_gene_anno_file = preorganized_snp_gene_annotation_dir + window_name + '_snp_gene_anno.npy'
		np.save(snp_gene_anno_file, variant_gene_anno)

		# Print to output
		t.write(line + '\t' + linked_genes_file + '\t' + snp_gene_anno_file + '\n')



	f.close()



t.close()


