'''The script shorten both 5' and 3' UTRs from gene models in BED format.'''

import sys
import csv

LEFT_UTR_LEN = 50
RIGHT_UTR_LEN = 50

class ExonObj(object):
    def __init__(self, chrom, start, end, trans_id, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = (end - start) + 1
        self.remove = False
        self.trans_id = trans_id
        self.strand = strand


def parse_bed(bed_file):
    '''Reads alignments from BED format and creates
    exon objects from a transcript.

    '''

    reader = csv.reader(open(bed_file), dialect='excel-tab')
    for row in reader:
        exons = []
        chrom = row[0]
        chrom_start = int(row[1])
        trans_id = row[3]
        cds_start = int(row[6])
        cds_end = int(row[7])
        strand = row[5]

        exon_sizes = [int(s) for s in row[10].split(',')]
        exon_starts = [chrom_start + int(s) for s in row[11].split(',')]

        for i in range(len(exon_starts)):
            exon_start = exon_starts[i]
            exon_end = exon_start + exon_sizes[i]

            exon = ExonObj(chrom, exon_start, exon_end, trans_id, strand)
            exons.append(exon)

        yield exons, cds_start, cds_end


def shorten_UTRs(exons, cds_start, cds_end):
    if (cds_start == exons[0].start and  # if no UTRs, do nothing
            cds_end == exons[-1].end):
        return

    cds_start_exon = None
    cds_end_exon = None
    for i in range(len(exons)):
        if (exons[i].start < cds_start) and (exons[i].end > cds_start):
            cds_start_exon = i
        if (exons[i].start < cds_end) and (exons[i].end > cds_end):
            cds_end_exon = i

    if cds_start_exon != None:
        j = cds_start_exon
        left_utr_length = cds_start - exons[j].start
        if left_utr_length > LEFT_UTR_LEN:
            exons[j].start = cds_start - LEFT_UTR_LEN
            # print 'shorten exon', exons[j].start, exons[j].end
        else:
            r = LEFT_UTR_LEN - left_utr_length
            while (left_utr_length < LEFT_UTR_LEN):
                if j > 0: j -= 1
                else: break

                if exons[j].end - r < exons[j].start:
                    left_utr_length += exons[j].end - exons[j].start
                    r = LEFT_UTR_LEN - left_utr_length
                else:
                    exons[j].start = exons[j].end - r
                    left_utr_length += r

        for i in range(0,j):
            exons[i].remove = True
            # print '***', exons[i].start, exons[i].end, exons[i].remove

    if cds_end_exon != None:
        # print 'j is ', j
        j = cds_end_exon
        right_utr_length = exons[j].end - cds_end

        # print 'right utr length', right_utr_length

        if right_utr_length > RIGHT_UTR_LEN:
            exons[j].end = cds_end + RIGHT_UTR_LEN
            # print 'new end', exons[j].end
        else:
            r = RIGHT_UTR_LEN - right_utr_length
            while (right_utr_length < RIGHT_UTR_LEN):
                if j < (len(exons) - 1): j += 1
                else: break

                if exons[j].start + r > exons[j].end:
                    right_utr_length += exons[j].end - exons[j].start
                    r = RIGHT_UTR_LEN - right_utr_length
                else:
                    exons[j].end = exons[j].start + r
                    right_utr_length += r
        for i in range(j + 1, len(exons)):
            exons[i].remove = True
            # print '***', exons[i].start, exons[i].end, exons[i].remove

    return


def print_bed(exons, cds_start, cds_end):
    '''Print a splice graph in BED format.'''

    chrom_start = exons[0].start
    chrom_end = exons[-1].end
    chrom = exons[0].chrom

    block_starts = ','.join([str(exon.start - chrom_start) for exon in exons])
    block_sizes = ','.join([str(exon.end - exon.start) for exon in exons])

    name = exons[0].trans_id
    strand = exons[0].strand
    score = 1000
    item_RGB = '0,0,0'
    block_count = len(exons)

    writer = csv.writer(sys.stdout, dialect='excel-tab')
    writer.writerow((chrom,
                    chrom_start,
                    chrom_end,
                    name,
                    score,
                    strand,
                    cds_start,
                    cds_end,
                    item_RGB,
                    block_count,
                    block_sizes,
                    block_starts))


def main():
    bedfile = sys.argv[1]
    for exons, cds_start, cds_end in parse_bed(bedfile):
        shorten_UTRs(exons, cds_start, cds_end)
        new_exons = []
        for exon in exons:
            if not exon.remove: new_exons.append(exon)

        print_bed(new_exons, cds_start, cds_end)


if __name__=='__main__':
    main()
