import sys

from pygr import seqdb

def getSeq(filename, genome):
    for line in open(filename):
        chrom, coord = line.strip().split(':')
        start, end = [int(i) for i in coord.split('-')]
        seq = genome[chrom][start-2:start+3]
        if 'N' not in seq:
            print str(seq).upper()


if __name__=='__main__':
    filename = sys.argv[1]
    genomeFile = sys.argv[2]
    genome = seqdb.SequenceFileDB(genomeFile, verbose=False)
    getSeq(filename, genome)
