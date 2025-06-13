import numpy as np
import os
import sys
import pdb
import gzip


def create_mapping_from_chrom_num_to_min_max_variant_position(sumstats_file):
	# Initialize dictionary
	dicti = {}
	for chrom_num in range(1,23):
		chrom_string = str(chrom_num)
		dicti[chrom_string] = (1000000000000, -10)

	for chrom_num in range(1,23):
		snp_file = baselineLD_anno_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'
		head_count = 0
		f = gzip.open(snp_file)
		for line in f:
			line = line.decode('utf-8').rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			# extract relevent fields
			chrom_string = data[0]
			snp_pos = int(data[1])
			old_tuple = dicti[chrom_string]
			dicti[chrom_string] = (min(old_tuple[0], snp_pos), max(old_tuple[1], snp_pos))
		f.close()
	return dicti


def window_does_not_overlap_long_ld_region(long_ld_regions, chrom_string, window_start, window_end):
	boolean = True
	for long_ld_region in long_ld_regions:
		long_ld_region_info = long_ld_region.split('_')
		long_ld_region_chrom = long_ld_region_info[0]
		long_ld_region_start = int(long_ld_region_info[1])
		long_ld_region_end = int(long_ld_region_info[2])
		if 'chr' + chrom_string != long_ld_region_chrom:
			continue

		# window start in the long ld region
		if window_start >= long_ld_region_start and window_start <= long_ld_region_end:
			boolean = False
		# window end in the long ld region
		if window_end >= long_ld_region_start and window_end <= long_ld_region_end:
			boolean = False
		if long_ld_region_start >= window_start and long_ld_region_start <= window_end:
			boolean = False
		if long_ld_region_end >= window_start and long_ld_region_end <= window_end:
			boolean = False

	return boolean





######################
# Command line args
######################
baselineLD_anno_dir = sys.argv[1]
quasi_independent_ld_blocks_file = sys.argv[2]
genome_wide_ld_windows_file = sys.argv[3]



# Regions identified in methods of Weissbrod et al 2020 Nature Genet (hg19) and lifted over to hg38
# We will throw out windows in these regions
long_ld_regions = ['chr6_25499772_33532223', 'chr8_8142478_12142491', 'chr11_45978449_57232526']  # hg38



# Print windows to output file
t = open(genome_wide_ld_windows_file, 'w')
t.write('chrom_num\tstart_pos_inclusive\tend_position_exclusive\n')


f = open(quasi_independent_ld_blocks_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	chrom_string = data[0].split('hr')[1]
	window_start = int(data[1])
	window_end = int(data[2])

	if window_does_not_overlap_long_ld_region(long_ld_regions, chrom_string, window_start, window_end):
		t.write(chrom_string + '\t' + str(window_start) + '\t' + str(window_end) + '\n')

f.close()
t.close()


