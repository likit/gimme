'''The script reads a BAM file and report SNPs
according to a position on a read.

Note for a coordinate in SAM/BAM format:
    In SAM file, it is 1-based.
    In BAM file, it is 0-based.
    (http://genome.sph.umich.edu/wiki/SAM)
'''

import sys
import pysam
# import matplotlib.pyplot as plot


class Read(object):
    def __init__(self, seq, tags):
        self.seq = seq
        self.tags = tags


def find_snp(aligned_read):
    '''Note for Cigar tag:
            I means a deletion in a read.
            D means an insertion in a read.
    '''

    mm = {}
    pos = 0
    n = ''
    INSERT = False

    md = dict(aligned_read.tags)['MD']

    for i in md:
        if i.isdigit():
            if INSERT:
                INSERT = False
            n += i
        else:
            if n:
                pos += int(n)
                n = ''
            if i.isalpha():
                if not INSERT:
                    pos += 1
                    try:
                        mm[pos] = aligned_read.seq[pos-1]
                    except IndexError:
                        print >> sys.stderr, aligned_read.seq, pos, md, mm
                        raise SystemExit
            else:
                if i == '^':
                    INSERT = True
    return mm


def read_bam(bamfile, reference=None):
    max_rlen = 0
    all_mismatches = {
                        'A':{},
                        'C':{},
                        'G':{},
                        'T':{},
                        'N':{},
                        '*':{}, # deletion
                        }

    print >> sys.stderr, 'Reading BAM file...'
    for n, aligned_read in enumerate(bamfile.fetch(
                                    reference), start=1):
        mismatches = find_snp(aligned_read)
        if aligned_read.rlen > max_rlen:
            max_rlen = aligned_read.rlen

        for pos, base in mismatches.iteritems():
            if pos not in all_mismatches[base]:
                all_mismatches[base][pos] = 1
            else:
                all_mismatches[base][pos] += 1

        if n % 1000 == 0:
            print >> sys.stderr , '...', n

    return all_mismatches, max_rlen


def plot_chart(all_mismatches, max_rlen):
    A_count = [all_mismatches['A'].get(i, 0) for i in range(1, max_rlen+1)]
    C_count = [all_mismatches['C'].get(i, 0) for i in range(1, max_rlen+1)]
    G_count = [all_mismatches['G'].get(i, 0) for i in range(1, max_rlen+1)]
    T_count = [all_mismatches['T'].get(i, 0) for i in range(1, max_rlen+1)]
    # N_count = [all_mismatches['N'].get(i, 0) for i in range(max_rlen)]
    # D_count = [all_mismatches['*'].get(i, 0) for i in range(max_rlen)]

    plot.plot(range(1, max_rlen+1), A_count, label='A')
    plot.plot(range(1, max_rlen+1), C_count, label='C')
    plot.plot(range(1, max_rlen+1), G_count, label='G')
    plot.plot(range(1, max_rlen+1), T_count, label='T')
    # plot.plot(range(max_rlen), N_count, label='N')
    # plot.plot(range(max_rlen), D_count, label='deletion')

    plot.xlabel('position')
    plot.ylabel('mismatch')
    plot.yscale('log')
    plot.legend()
    plot.show()


def write_output(all_mismatches, max_rlen):
    A_count = [str(all_mismatches['A'].get(i, 0)) for i in range(1, max_rlen+1)]
    C_count = [str(all_mismatches['C'].get(i, 0)) for i in range(1, max_rlen+1)]
    G_count = [str(all_mismatches['G'].get(i, 0)) for i in range(1, max_rlen+1)]
    T_count = [str(all_mismatches['T'].get(i, 0)) for i in range(1, max_rlen+1)]

    print '# Rows are A,C,G,T.'
    print ','.join(A_count)
    print ','.join(C_count)
    print ','.join(G_count)
    print ','.join(T_count)


def main():
    bamfile = pysam.Samfile(sys.argv[1], 'rb')

    try:
        reference = sys.argv[2]
    except IndexError:
        reference = None

    # test_find_snp()

    mismatches, max_rlen = read_bam(bamfile, reference)
    # plot_chart(mismatches, max_rlen)
    write_output(mismatches, max_rlen)


def test_find_snp():
    '''Test one insertion.'''
    seq = "ACAAAAAGGAGACCTCTTTCTTCAGGAAAAAAAAAAAGCCTTCATTTCCCCTTCATCTCTTTGTGCTGCCATAAC"
    tags = [('MD', '3G22^A48A0')]

    read = Read(seq, tags)

    mm = find_snp(read)
    print mm

    # tags = '44^A21'
    # mm = find_snp(tags)
    # assert mm == {45:'*'}, "Failed %s != {45:'*'}" % mm

    # '''Test two insertions.'''
    # tags = '7^TG58'
    # mm = find_snp(tags)
    # assert mm == {8:'*', 9:'*'}, "Failed %s != {8:'*', 9:'*'}" % mm

    # '''Test ending polymorphism.'''
    # tags = '64G0'
    # mm = find_snp(tags)
    # assert mm == {65:'G'}, "Failed %s != {65:'G'}" % mm

    # '''Test one insertion + snp.'''
    # tags = '54^C3A0'
    # mm = find_snp(tags)
    # assert mm == {55:'*', 59:'A'}, "Failed %s != {55:'*', 59:'A'}" % mm

    # '''Test adjacent insertion and snp.'''
    # tags = '61^C0G3'
    # mm = find_snp(tags)
    # assert mm == {62:'*', 63:'G'}, "Failed %s != {62:'*', 63:'G'}" % mm


if __name__=='__main__':
    main()
