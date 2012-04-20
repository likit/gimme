'''The script reads a sequence in FASTA format and
identify a repeat region due to misassembly (a bubble)
and collapse it.

Author: Likit Preeyanon
Email: preeyano@msu.edu
'''

import sys
from collections import namedtuple

KMERSIZE = 27

Sequence = namedtuple('Sequence', ['id', 'seq'])

def parse_fasta(filename):
    seqid = None
    for line in open(filename):
        if line.startswith('>'):
            newid = line.strip().lstrip('>')
            if seqid == None:
                seqid = newid
                seq = ''
            elif newid != seqid:
                yield Sequence(seqid, seq)
                seqid = newid
                seq = ''
        else:
            seq += line.strip()

    yield Sequence(seqid, seq)


def makehash(sequence, kmersize=KMERSIZE):
    kmer_table = {}
    all_kmers = []
    if sequence.seq < kmersize:
        return None
    
    kmerid = 0
    end = kmersize
    kmer = sequence.seq[:end]
    kmer_table[kmer] = [kmerid]
    all_kmers.append(kmer)
    start = 1
    end += 1
    while len(sequence.seq) - start >= kmersize:
        kmerid += 1
        kmer = sequence.seq[start:end]
        try:
            kmer_table[kmer].append(kmerid)
        except KeyError:
            kmer_table[kmer] = [kmerid]

        all_kmers.append(kmer)

        #print >> sys.stderr, kmer, kmer_table[kmer]

        start += 1
        end += 1

    the_rest = sequence.seq[start:]

    return kmer_table, all_kmers, the_rest


def collapse(kmer_table, all_kmers):
    discard = set()
    for kmer in all_kmers:
        if len(kmer_table[kmer]) > 0:
            for node in kmer_table[kmer][1:]:
                discard.add(node - 1)
                discard.add(node)

    return discard


def rebuild_sequence(kmer_table, all_kmers, discard):
    sequence = all_kmers[0]
    for kmer in all_kmers[1:]:
        if not set(kmer_table[kmer]).intersection(discard):
            sequence += kmer[-1]
        else:
            continue

    return sequence

def main(argv):
    fasta_file = argv[1]

    try:
        kmersize = int(argv[2])
    except IndexError:
        kmersize = None

    log_file = open(fasta_file+'.log', 'w')

    for n, sequence in enumerate(parse_fasta(fasta_file), start=1):
        kmer_table, all_kmers, the_rest = makehash(sequence, kmersize)
        discard = collapse(kmer_table, all_kmers)
        new_sequence = rebuild_sequence(kmer_table, all_kmers, discard)
        print '>%s\n%s' % (sequence.id, new_sequence)
        if discard:
            print >> log_file, '%s\t%d\t%d'\
                    % (sequence.id, len(sequence.seq), len(discard))

        if n % 1000 == 0:
            print >> '...', n

    log_file.close()


if __name__=='__main__':
    main(sys.argv)
