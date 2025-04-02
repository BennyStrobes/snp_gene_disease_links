import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np
import os
import pdb
import gzip









##################
# Command line args
###################
ld_matrix_output_dir = sys.argv[1]

# Create new aggregate output file
output_file = ld_matrix_output_dir + 'window_LD_summary_cross_chromosomes.txt'
t = open(output_file,'w')

for chrom_num in range(1,23):
	filer = ld_matrix_output_dir + 'window_LD_summary_chr' + str(chrom_num) + '.txt'
	f = open(filer)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			if chrom_num == 1:  # only print to output on first chromosome
				t.write(line + '\n')
			continue
		t.write(line + '\n')
	f.close()
t.close()
