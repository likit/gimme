'''The script gives a number of genes from BED file.
A gene name must be in "chromosome:gene.isoform" format.
'''

import sys, csv

genes = set([])

for rows in csv.reader(open(sys.argv[1]), dialect='excel-tab'):
    gene, isoform = rows[3].split('.')
    genes.add(gene)

print 'total genes = %d' % len(genes)
