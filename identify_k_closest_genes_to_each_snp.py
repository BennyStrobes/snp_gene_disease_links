import numpy as np
import os
import sys
import pdb
import gzip





def extract_snps_on_this_chromosome(baselineLD_anno_dir, chrom_num):
	snp_file = baselineLD_anno_dir + 'baselineLD.' + chrom_num + '.annot.gz'
	head_count = 0
	rsids = []
	pos = []
	f = gzip.open(snp_file)
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# extract relevent fields
		line_chrom_num = data[0]
		snp_pos = int(data[1])
		rsid = data[2]
		rsids.append(rsid)
		pos.append(snp_pos)
		if line_chrom_num != chrom_num:
			print('assumptioner ororr')
			pdb.set_trace()

	f.close()

	rsids = np.asarray(rsids)
	pos = np.asarray(pos)
	if len(pos) != len(rsids):
		print('assumption eororro')
		pdb.set_trace()
	return rsids, pos


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


######################
# Command line args
######################
baselineLD_anno_dir = sys.argv[1]
gene_annotation_file = sys.argv[2]
K_closest_genes = int(sys.argv[3])
chrom_num = sys.argv[4]
snp_gene_annotation_dir = sys.argv[5]



#############
# First extract snps on this chromosome
snp_rsids, snp_positions = extract_snps_on_this_chromosome(baselineLD_anno_dir, chrom_num)

#############
# Second extract snps on this chromosome
gene_ensamble_ids, gene_positions = create_mapping_from_gene_name_to_gene_info(gene_annotation_file, chrom_num)


output_file = snp_gene_annotation_dir + 'snp_gene_pairs_chr' + str(chrom_num) + '.txt'
t = open(output_file,'w')
t.write('snp_id\tsnp_chrom\tsnp_pos\tgenes\tgene_distances\n')


valid_genes = {}
# loop through snps
for ii, snp_rsid in enumerate(snp_rsids):
	snp_position = snp_positions[ii]

	# Calculate distance between this snp and all genes
	distances = np.abs(snp_position - gene_positions)

	# Ordered closest gene indices
	ordered_closest_gene_indices = np.argsort(distances)[:K_closest_genes]

	# Extract names and distances of K closest genes
	ordered_closest_gene_names = gene_ensamble_ids[ordered_closest_gene_indices]
	ordered_closest_gene_distances = distances[ordered_closest_gene_indices]

	for gene in ordered_closest_gene_names:
		valid_genes[gene] = 1

	t.write(snp_rsid + '\t' + chrom_num + '\t' + str(snp_position) + '\t' + ';'.join(ordered_closest_gene_names) + '\t' + ';'.join(ordered_closest_gene_distances.astype(str)) + '\n')
t.close()


output_file = snp_gene_annotation_dir + 'genes_chr' + str(chrom_num) + '.txt'
t = open(output_file,'w')
t.write('gene_id\tchrom_num\ttss\n')
for ii, ensamble_id in enumerate(gene_ensamble_ids):
	if ensamble_id not in valid_genes:
		continue
	t.write(ensamble_id + '\t' + chrom_num + '\t' + str(gene_positions[ii]) + '\n')
t.close()


print(len(valid_genes))
print(len(gene_ensamble_ids))

