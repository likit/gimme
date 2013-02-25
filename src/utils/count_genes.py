'''The script gives a number of genes from BED file.
A gene name must be in "chromosome:gene.isoform" format.
'''

import sys
import csv

genes = set([])
isoforms = set([])

for rows in csv.reader(open(sys.argv[1]), dialect='excel-tab'):
    try:
        gene, isoform = rows[3].split('.')
        genes.add(gene)
        isoforms.add(rows[3])
    except:
        pass

print 'Total genes = %d' % len(genes)
print 'Total isoforms = %d' % len(isoforms)
