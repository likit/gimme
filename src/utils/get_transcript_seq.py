'''This script reads a gene model from BED file
and writes a DNA sequence to standard output.
The script requires pygr package.

'''

import sys
import csv

from collections import namedtuple
from pygr import seqdb, sequtil

Exon = namedtuple('Exon', 'chrom, start, end')

def get_sequence_transcript(genome, exons, strand='positive'):
    seq = ''
    for exon in exons:
        s = genome[exon.chrom][exon.start:exon.end]
        seq += str(s)
    if strand == 'positive':
        return seq
    elif strand == 'negative':
        return s.reverse_complement(seq)
    else:
        raise ValueError

def get_sequence_exon(genome, exons, strand='positive'):
    seqs = []
    for exon in exons:
        s = genome[exon.chrom][exon.start:exon.end]
        seqs.append(str(s))

    if strand == 'positive':
        return seqs
    elif strand == 'negative':
        rev_seqs = [s.reverse_complement(seq) for seq in seqs]
        return rev_seqs
    else:
        raise ValueError

def write_seq(filename, genome, output, strand):
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
        strand = 'negative' if line[5] == '-' else 'positive'

        if output == 'transcript':
            seq = get_sequence_transcript(genome, exons, strand)
            sequtil.write_fasta(sys.stdout, seq, id=gene_id)
        elif output == 'exon':
            seqs = get_sequence_exon(genome, exons, strand)

            for n, seq in enumerate(seqs, start=1):
                seq_id = gene_id + '_' + str(n)
                sequtil.write_fasta(sys.stdout, seq, id=seq_id)
        else:
            print >> sys.stderr, 'Unsupported output format.'
            raise SystemExit

        if n % 1000 == 0: print >> sys.stderr, '...', n

if __name__=='__main__':
    if (len(sys.argv) == 1 or sys.argv[1] == '-h'):
        print >> sys.stderr, \
            'Usage: python get_transcript_seq.py <bed file> ' + \
            '<genome file> [option]'
        print >> sys.stderr, \
            'Output options:\n' + \
            '\t-t\toutput a transcript sequence per record [default]\n' + \
            '\t-e\toutput an exon sequence per record\n' + \
            '\t-r\toutput reverse complement.\n' + \
            '\t-h\tprint help message.'
        raise SystemExit

    filename = sys.argv[1]
    genome_file = sys.argv[2]
    strand = 'positive'
    output = 'transcript'
    if len(sys.argv) > 3: # options specified
        for opt in sys.argv[3:]:
            if opt == '-e':
                output = 'exon'
            elif opt == '-t':
                pass
            elif opt == '-r':
                strand = 'negative'
            else:
                print >> sys.stderr, 'Unrecognized option %s' % opt
                raise SystemExit

    # print >> sys.stderr, filename, genome_file, output, strand
    genome = seqdb.SequenceFileDB(genome_file, verbose=False)
    write_seq(filename, genome, output, strand)
