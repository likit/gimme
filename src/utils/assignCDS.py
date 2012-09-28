'''The script reads gene models in BED format and translated proteins
in FASTA from ESTScan and reports gene models in BED format with ORFs
according to corresponding translated proteins.

'''
import sys, csv
from bx.intervals.intersection import Interval, Intersecter


class ExonObj(object):
    def __init__(self, seqId, start, end, cStart, cEnd, size):
        self.seqId = seqId
        self.start = start
        self.end = end
        self.cStart = cStart
        self.cEnd = cEnd
        self.geneSize = size

    def __str__(self):
        return '%s:%d-%d' % (self.seqId, self.start, self.end)


def parseBed(filename):
    '''Reads BED file and returns a group of exons
    of a transcript.

    '''

    with open(filename) as fp:
        for row in csv.reader(fp, dialect='excel-tab'):
            exons = []
            seqId = row[3]
            chromStart = int(row[1])

            blockSizes = [int(i) for i in row[-2].split(',')]
            chromStarts = [int(i) + chromStart
                                for i in row[-1].split(',')]
            chromEnds = [int(blockSizes[i]) + chromStarts[i]\
                                for i in range(len(blockSizes))]
            size = sum(blockSizes)

            i = 0
            start = i
            while i in range(len(blockSizes)):
                end = start + blockSizes[i]
                exons.append(ExonObj(seqId,
                                        start,
                                        end,
                                        chromStarts[i],
                                        chromEnds[i],
                                        size),
                                        )
                i += 1
                start = end

            yield row, exons


def parseCDS(filename):
    cds = {}
    for line in open(filename):
        if line.startswith('>'):
            line = line.strip('>')
            strand = '-' if 'minus' in line else '+'
            values = line.split()
            name = values[0].strip('>').strip(';')

            start = int(values[2]) - 1
            end = int(values[3]) - 1

            assert start < end, line.strip()
            cds[name] = (name, start, end, strand)

    return cds


def buildIntervalTree(exons):
    '''Build interval tree from exon annotations.'''
    tree = Intersecter()
    for exon in exons:
        tree.add_interval(Interval(exon.start, exon.end,
                            value={'cStart':exon.cStart,
                                    'cEnd':exon.cEnd}))
    return tree

def findCDS(bedFile, cds):
    """docstring for findCDS"""

    writer = csv.writer(sys.stdout, dialect='excel-tab')
    for row, exons in parseBed(bedFile):
        tree = buildIntervalTree(exons)

        trans_size = sum([(exon.cEnd - exon.cStart) for exon in exons])

        exon = exons[0]
        try:
            name, start, end, strand = cds[exon.seqId]
            if strand == '-':
                new_end = trans_size - start
                new_start = trans_size - end
                start = new_start
                end = new_end
                assert start < end, "%d-%d %d %s" % \
                            (start, end, trans_size, exon.seqId)

        except KeyError:
            row[5] = '.'
            writer.writerow(row)
        else:
            try:
                exon1 = tree.find(start, start + 1)[0]
                exon2 = tree.find(end - 1, end)[0]
            except IndexError:
                print >> sys.stderr, name, start, end, strand
            else:
                CDSStart = exon1.value['cStart'] + (start - exon1.start)
                CDSEnd = exon2.value['cStart'] + (end - exon2.start)
                row[5] = strand
                row[6] = CDSStart
                row[7] = CDSEnd
                writer.writerow(row)


if __name__=='__main__':
    print >> sys.stderr, 'Parsing FASTA file...'
    cds = parseCDS(sys.argv[1])
    findCDS(sys.argv[2], cds)
