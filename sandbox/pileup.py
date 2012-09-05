import sys
import csv
from collections import namedtuple

from pysam import Samfile

Exon = namedtuple('Exon', 'name, chrom, start, end')
uniqSeqs = set([])

class Exon(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = self.__str__()
        self.transcripts = []
        self.sequence = ''

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)

def getSequence(genome, exons):
    seq = ''
    for exon in exons:
        s = genome[exon.chrom][exon.start:exon.end]
        seq += str(s)

    return seq

def parseBed(bedFile):

    reader = csv.reader(open(filename), dialect='excel-tab')

    for n, line in enumerate(reader, start=1):
        chrom = line[0]
        chromStart = int(line[1])
        geneId = line[3]
        exonStarts = [int(start) + chromStart for start in line[-1].split(',')]
        exonSizes = [int(size) for size in line[-2].split(',')]
        exonEnds = [exonStarts[i] + exonSizes[i] for i in range(len(exonStarts))]
        exons = [Exon(chrom, exonStarts[i], exonEnds[i]) for i in \
                range(len(exonStarts))]

        #seq = getSequence(genome, exons)
        #sequtil.write_fasta(sys.stdout, seq, id=geneId)

        if n % 1000 == 0: print >> sys.stderr, '...', n

if __name__=='__main__':
    bedFile = sys.argv[1]
    #genome = seqdb.SequenceFileDB(sys.argv[2], verbose=False)
    parseBed(bedFile)
