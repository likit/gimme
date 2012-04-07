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
    try:
        for exon in exons:
            s = genome[exon.chrom][exon.start:exon.end]
            seq += str(s)
    except IndexError as e:
        print >> sys.stderr, exons
        raise e

    return seq

def parse_seq(filename, genome):
    reader = csv.reader(open(filename), dialect='excel-tab')

    for line in reader:
        chrom = line[13]
        gene_id = line[9]
        exon_starts = [int(start) for
                        start in line[-1].split(',')[:-1]]

        exon_sizes = [int(size) for size in line[-3].split(',')[:-1]]

        exon_ends = [exon_starts[i] + exon_sizes[i]
                        for i in range(len(exon_starts))]

        exons = [Exon(
                        chrom, exon_starts[i],
                        exon_ends[i])
                        for i in range(len(exon_starts))
                        ]
        print >> sys.stderr, exons, gene_id
        yield exons, gene_id


def main():
    filename = sys.argv[1]
    genome = seqdb.SequenceFileDB(sys.argv[2], verbose=False)
    for n, (exons, gene_id) in enumerate(
                    parse_seq(filename, genome), start=1):

        seq = get_sequence(genome, exons)
        sequtil.write_fasta(sys.stdout, seq, id=gene_id)

        if n % 1000 == 0: print >> sys.stderr, '...', n

if __name__=='__main__':
    main()
