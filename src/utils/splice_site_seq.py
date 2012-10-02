'''The script reads in splice junctions and retrieve sequences from
a given genome.

Splice junctions are obtained from compare_junctions.py with -a option.

Then run this command to filter out introns with ambiguous strand.

cat junctions.txt | sort -k 1 | uniq -c | awk '$1==1' > junctions.nr.txt

Junctions without a transcriptional direction are ignored.

'''

import sys
from pygr import seqdb


def parse_input(infile):
    with open(infile) as fp:
        for line in fp:
            try:
                _, junction, strand = line.split()
            except ValueError:
                junction, strand = line.split()

            chrom, coord = junction.split(':')
            start, end = coord.split('-')
            start = int(start)
            end = int(end)

            if strand == '.': continue

            yield chrom, start, end, strand


def get_sequence(genome, intron):
    chrom, start, end, strand = intron
    ref = genome[chrom]
    if strand == '+':
        donor_start = start - 5 if start - 5 > 0 else None
        donor_end = start + 5 if start + 5 < len(ref) else None
        if donor_start and donor_end:
            donor_seq = ref[donor_start:donor_end]

        acceptor_start = end - 5 if end - 5 > 0 else None
        acceptor_end = end + 5 if end + 5 < len(ref) else None
        if acceptor_start and acceptor_end:
            acceptor_seq = ref[acceptor_start:acceptor_end]

        return str(donor_seq), str(acceptor_seq)

    if strand == '-':
        acceptor_start = start - 5 if start - 5 > 0 else None
        acceptor_end = start + 5 if start + 5 < len(ref) else None
        if acceptor_start and acceptor_end:
            acceptor_seq = ref[acceptor_start:acceptor_end]
            acceptor_seq = acceptor_seq.reverse_complement(str(acceptor_seq))

        donor_start = end - 5 if end - 5 > 0 else None
        donor_end = end + 5 if end + 5 < len(ref) else None
        if donor_start and donor_end:
            donor_seq = ref[donor_start:donor_end]
            donor_seq = donor_seq.reverse_complement(str(donor_seq))

        return str(donor_seq), str(acceptor_seq)


def main():
    infile = sys.argv[1]
    refseq = sys.argv[2]

    refseq = seqdb.SequenceFileDB(refseq)
    op1 = open('donor_sites', 'w')
    op2 = open('acceptor_sites', 'w')
    for n, intron in enumerate(parse_input(infile), start=1):
        intron_str = '%s:%d-%d' % intron[:-1]
        donor, acceptor = get_sequence(refseq, intron)
        print >> op1, '>%s\n%s' % (intron_str, donor)
        print >> op2, '>%s\n%s' % (intron_str, acceptor)

        if n % 1000 == 0: print >> sys.stderr, '...', n

    op1.close()
    op2.close()


if __name__ == '__main__':
    main()
