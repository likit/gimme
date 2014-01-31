'''This script converts GTF/GFF format to BED format.
The output is written to standard output.

Usage gff2bed.py file.gtf/gff

'''

import csv
import sys
from collections import namedtuple

stderr = sys.stderr

Transcript = namedtuple('Transcript',
        ['chrom', 'id', 'strand', 'exons', 'geneid'])


def parse(filename):
    reader = csv.reader(open(filename), dialect='excel-tab')
    exons = []
    name = None
    strand = None
    for row in reader:
        if row[2] == 'exon':
            transId = eval(row[-1].split('; ')[1].split(' ')[-1])
            geneId = eval(row[-1].split('; ')[0].split(' ')[-1])
            #print transId, geneId

            if not name:
                name = transId
                gene_name = geneId

            if transId == name:
                chrom = row[0]
                start = int(row[3]) - 1
                end = int(row[4])
                strand = row[6]
                exons.append((start, end))
            else:
                exons = sorted(exons, key=lambda x: x[0])
                yield Transcript(chrom, name, strand, exons, gene_name)

                chrom = row[0]
                start = int(row[3]) - 1
                end = int(row[4])
                strand = row[6]
                name = transId
                gene_name = geneId
                exons = [(start, end)]

    yield Transcript(chrom, name, strand, exons, gene_name)


def printBED(transcript):
    blockStarts = []
    blockSizes = []
    chromStart = transcript.exons[0][0]

    for start, end in transcript.exons:
        blockStarts.append(start - chromStart)
        blockSizes.append(end - start)

    chromEnd = end

    writer = csv.writer(sys.stdout, dialect='excel-tab')
    writer.writerow([transcript.chrom,
                    chromStart,
                    chromEnd,
                    "%s.%s" % (transcript.geneid, transcript.id),
                    1000,
                    transcript.strand,
                    chromStart,
                    chromEnd,
                    '255,0,0',
                    len(blockStarts),
                    ','.join([str(s) for s in blockSizes]),
                    ','.join([str(s) for s in blockStarts])])


if __name__ == '__main__':
    for transcript in parse(sys.argv[1]):
        printBED(transcript)
