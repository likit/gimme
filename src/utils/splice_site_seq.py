'''The script reads in splice junctions and retrieve sequences from
a given genome.

Splice junctions are obtained from compare_junctions.py with -a option.
Then run this command to filter out introns with ambiguous strand.

cat junctions.txt | sort -k 1 | uniq -c | awk '$1==1' > junctions.nr.txt

Junctions without a transcriptional direction are ignored.

'''

import sys
from pygr import seqdb

SEQLEN = 21  # number of nucleotides on each side of the splice junction.


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

            if strand == '.':
                continue

            yield chrom, start, end, strand


def get_splice_site(genome, intron):
    chrom, start, end, strand = intron
    ref = genome[chrom]
    if strand == '+':
        donor_start = start - SEQLEN if start - SEQLEN > 0 else None
        donor_end = start + SEQLEN if start + SEQLEN < len(ref) else None
        if donor_start and donor_end:
            donor_seq = ref[donor_start:donor_end]

        acceptor_start = end - SEQLEN if end - SEQLEN > 0 else None
        acceptor_end = end + SEQLEN if end + SEQLEN < len(ref) else None
        if acceptor_start and acceptor_end:
            acceptor_seq = ref[acceptor_start:acceptor_end]

        return str(donor_seq), str(acceptor_seq)

    if strand == '-':
        acceptor_start = start - SEQLEN if start - SEQLEN > 0 else None
        acceptor_end = start + SEQLEN if start + SEQLEN < len(ref) else None
        if acceptor_start and acceptor_end:
            acceptor_seq = ref[acceptor_start:acceptor_end]
            acceptor_seq = acceptor_seq.reverse_complement(str(acceptor_seq))

        donor_start = end - SEQLEN if end - SEQLEN > 0 else None
        donor_end = end + SEQLEN if end + SEQLEN < len(ref) else None
        if donor_start and donor_end:
            donor_seq = ref[donor_start:donor_end]
            donor_seq = donor_seq.reverse_complement(str(donor_seq))

        return str(donor_seq), str(acceptor_seq)


def get_junction_seq(genome, intron):
    chrom, start, end, strand = intron
    ref = genome[chrom]

    donor_end = start - 1
    donor_start = donor_end - SEQLEN

    if donor_start and donor_end:
        donor_seq = ref[donor_start:donor_end]

    acceptor_start = end + 1
    acceptor_end = acceptor_start + SEQLEN

    if acceptor_start and acceptor_end:
        acceptor_seq = ref[acceptor_start:acceptor_end]

    if strand == '+':
        return str(donor_seq) + str(acceptor_seq)
    else:
        donor_seq = donor_seq.reverse_complement(str(donor_seq))
        acceptor_seq = acceptor_seq.reverse_complement(str(acceptor_seq))
        return str(acceptor_seq) + str(donor_seq)


def main():
    infile = sys.argv[1]
    refseq = sys.argv[2]

    refseq = seqdb.SequenceFileDB(refseq)
    # op1 = open('donor_sites', 'w')
    # op2 = open('acceptor_sites', 'w')
    for n, intron in enumerate(parse_input(infile), start=1):
        intron_str = '%s:%d-%d' % intron[:-1]
        # donor, acceptor = get_splice_site(refseq, intron)
        # print >> op1, '>%s\n%s' % (intron_str, donor)
        # print >> op2, '>%s\n%s' % (intron_str, acceptor)
        junction_seq = get_junction_seq(refseq, intron)
        print '>%s\n%s' % (intron_str, junction_seq)

        if n % 1000 == 0:
            print >> sys.stderr, '...', n

    # op1.close()
    # op2.close()


if __name__ == '__main__':
    main()
