'''This script reads a gene model from PSL file
and writes a DNA sequence to standard output.
The script requires pygr package.

'''

import sys
import csv

from collections import namedtuple
from pygr import seqdb, sequtil

Exon = namedtuple('Exon', 'chrom, start, end')

def get_sequence(genome, exons):
    seq = ''
    for exon in exons:
        s = genome[exon.chrom][exon.start:exon.end]
        seq += str(s)

    return seq

def write_seq(filename, genome):
    reader = csv.reader(open(filename), dialect='excel-tab')

    for n, line in enumerate(reader, start=1):
        chrom = line[0]
        chrom_start = int(line[1])
        gene_id = line[3]
        exon_starts = [int(start) + chrom_start for start in line[-1].split(',')]
        exon_sizes = [int(size) for size in line[-2].split(',')]
        exon_ends = [exon_starts[i] + exon_sizes[i] for i in range(len(exon_starts))]
        exons = [Exon(chrom, exon_starts[i], exon_ends[i]) for i in \
                range(len(exon_starts))]

        seq = get_sequence(genome, exons)
        sequtil.write_fasta(sys.stdout, seq, id=gene_id)

        if n % 1000 == 0: print >> sys.stderr, '...', n

if __name__=='__main__':
    filename = sys.argv[1]
    genome = seqdb.SequenceFileDB(sys.argv[2], verbose=False)
    write_seq(filename, genome)
