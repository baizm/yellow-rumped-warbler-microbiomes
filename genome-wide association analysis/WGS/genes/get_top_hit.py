##--script to take first alignment from blastn output file, format---7##

import sys
import re

in_file = sys.argv[1]
out_file = sys.argv[2]

sys.stdout=open(out_file, 'w')

print('qseqid'+'\t'+'sseqid'+'\t'+'pident'+'\t'+'length'+'\t'+'mismatch'+'\t'+'gapopen'+'\t'+'qstart'+'\t'+'qend'+'\t'+'sstart'+'\t'+'send'+'\t'+'evalue'+'\t'+'bitscore')

with open(in_file, 'r') as f:
    for line in f:
        if 'hits found' in line:
        	split = re.split(' ', line)
        	n = split[1]
        	if int(n) > 0:
        		next_line = f.readline()
        		print(next_line,end='')

#usage: python get_top_hit.py in.txt out.txt

#ran on ccr:
#python get_top_hit.py taegut_genes_a.txt taegut_genes_top_a.txt
#python get_top_hit.py taegut_genes_b1.txt taegut_genes_top_b1.txt
#python get_top_hit.py taegut_genes_b2.txt taegut_genes_top_b2.txt
#python get_top_hit.py taegut_genes_all.txt taegut_genes_top_all.txt
