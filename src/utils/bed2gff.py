'''The script converts BED format to GFF format.
Output is written to standrad output.

'''

import sys
import csv


def parseBED(filename):
    reader = csv.reader(open(filename), dialect='excel-tab')
    for row in reader:
        chrom = row[0]
        chrom_start = int(row[1]) + 1
        transcript_id = row[3]
        strand = row[5]
        exon_sizes = [int(s) for s in row[10].split(',')]
        exon_starts = [int(s) for s in row[11].split(',')]
        yield (chrom,
                chrom_start,
                transcript_id,
                strand,
                exon_sizes,
                exon_starts)


def printGFF(transcript, source="custom"):
    chrom = transcript[0]
    chrom_start = transcript[1]
    transcript_id = transcript[2]
    gene_id = transcript_id.split(".")[0]
    strand = transcript[3]
    exon_sizes = transcript[4]
    exon_starts = transcript[5]

    exons = []

    for i in range(len(exon_sizes)):
        exon_start = exon_starts[i] + chrom_start
        exon_end = exon_start + exon_sizes[i] - 1
        exons.append((exon_start, exon_end))

    exon_id = 0
    for exon in exons:
        exon_id += 1
        attributes = 'gene_id \"%s\"; transcript_id \"%s\"; '
        attributes += 'exon_number \"%d\"; gene_name \"%s\"; '
        attributes += 'transcript_name \"%s\"'
        attributes = attributes % \
                (gene_id, transcript_id, exon_id, gene_id, transcript_id)
        print "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s" % (chrom,
                                                        source,
                                                        "exon",
                                                        exon[0],
                                                        exon[1],
                                                        ".",
                                                        strand,
                                                        ".",
                                                        attributes)


if __name__ == '__main__':
    filename = sys.argv[1]
    if len(sys.argv) == 2:
        source = "custom"
    else:
        source = sys.argv[2]
    for transcript in parseBED(filename):
        printGFF(transcript, source)
