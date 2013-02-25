'''The script converts alignments in PSL format to BED format. '''

import sys
import csv
import pslparser


def parsePSL(filename):
    for pslObj in pslparser.read(open(filename)):
        strand = pslObj.strand
        chrom = pslObj.tName
        name = pslObj.qName
        chromStart = pslObj.tStart
        chromEnd = pslObj.tEnd
        blockCount = pslObj.blockCount
        blockSizes = ','.join([str(i) for i in pslObj.blockSizes])
        blockStarts = ','.join([str(i - chromStart) \
                                        for i in pslObj.tStarts])
        thickStart = chromStart
        thickEnd = chromEnd

        writer = csv.writer(sys.stdout, dialect='excel-tab')
        writer.writerow((chrom,
                        chromStart,
                        chromEnd,
                        name,
                        255,
                        strand,
                        thickStart,
                        thickEnd,
                        '0,0,0',
                        blockCount,
                        blockSizes,
                        blockStarts))

if __name__ == '__main__':
    parsePSL(sys.argv[1])
