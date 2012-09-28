''' This script reports differential splice junctions 
of two gene models. It can read models in PSL and BED format.
Please adjust the parser function according to the data format.

All splice junctions from input files are written in first_sp.txt and
second_sp.txt respectively.

Differential introns are written in first_diff_sp.txt and
second_diff_sp.txt.

These files will not be created if there is no differential splice
junctions.

'''


import csv, os.path, sys
from collections import namedtuple

import pslparser

MIN_INTRON = 50
MAX_INTRON = 300000

Transcript = namedtuple('Transcript', [
                                        'chrom',
                                        'start',
                                        'end',
                                        'blockSizes',
                                        'blockStarts',
                                        'name',
                                        'strand',
                                    ])

class Exon(object):
    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.set_name()
        self.introns = set()
        self.terminal = None
        self.next_exons = set()
        self.prevExons = set()
    
    def set_name(self):
        self.name = self.__str__()

    def __str__(self):
        return '%s:%d-%d\t%s' % (self.chrom, self.start,
                                        self.end, self.strand)


def delete_gap(exons):
    i = 0

    new_exon = []
    curr_exon = exons[i]

    while True:
        try:
            next_exon = exons[i + 1]
        except IndexError:
            break
        else:
            if next_exon.start - curr_exon.end <= MIN_INTRON:
                curr_exon.end = next_exon.end
            else:
                new_exon.append(curr_exon)
                curr_exon = next_exon

            curr_exon.set_name()
        i += 1

    new_exon.append(curr_exon)
    return new_exon


def parse_psl(filename):
    '''Reads alignments from PSL format and create
    exon objects from each transcript.

    '''

    for pslobj in pslparser.read(open(filename)):
        exons = []
        for i in range(len(pslobj.tStarts)):
            start = pslobj.tStarts[i]
            end = start + pslobj.blockSizes[i]
            strand = pslobj.strand

            exon = Exon(pslobj.tName, start, end, strand)
            exons.append(exon)
        yield exons


def add_intron(exons, intron_db, all=False):
    for i in range(len(exons)):
        curr_exon = exons[i]
        try:
            next_exon = exons[i + 1]
        except IndexError:
            continue
        else:
            intron_start = curr_exon.end + 1
            intron_end = next_exon.start - 1
            intron_strand = curr_exon.strand
            intron_chrom = curr_exon.chrom
            intron_size = intron_end - intron_start
            if all:
                intron = "%s:%d-%d\t%s" % (intron_chrom,
                                                intron_start,
                                                intron_end,
                                                intron_strand,
                                            )
            else:
                intron = "%s:%d-%d" % (intron_chrom,
                                            intron_start,
                                            intron_end,
                                        )

            try:
                assert intron_start < intron_end
            except AssertionError:
                pass
            else:
                if (intron_size > MIN_INTRON
                        and intron_size < MAX_INTRON):
                    intron_db.add(intron)


def parse_bed(bedFile):
    with open(bedFile) as fp:
        for row in csv.reader(fp, dialect='excel-tab'):
            blockStarts = [int(i) for i in row[11].split(',')]
            blockSizes = [int(i) for i in row[10].split(',')]
            strand = row[5]
            transcript = Transcript(row[0],
                                        int(row[1]),
                                        int(row[2]),
                                        blockSizes,
                                        blockStarts,
                                        row[3],
                                        strand,
                                    )
            yield transcript


def add_intron_bed(transcript, modelintron_db, all=False):
    exons = []

    for i in range(len(transcript.blockStarts)):
        start = transcript.blockStarts[i] + transcript.start
        end = transcript.blockSizes[i] + start
        strand = transcript.strand
        exon = Exon(transcript.chrom, start, end, strand)

        exons.append(exon)

    add_intron(exons, modelintron_db, all)


def adjust_parser(inputFile1, format1,
                    inputFile2, format2, db1, db2):

    parsers = {'psl':parse_psl, 'bed':parse_bed}

    parser1 = parsers[format1]
    parser2 = parsers[format2]
    first = (inputFile1, parser1, db1)
    second = (inputFile2, parser2, db2)

    return first, second

def main(args):
    db1 = set()
    if args.all:  # Print out all splice junctions.
        if args.psl:
            parser = parse_psl
            input_filename1 = args.psl[0]
        elif args.bed:
            parser = parse_bed
            input_filename1 = args.bed[0]
        else:  # default format
            parser = parse_bed
            input_filename1 = args.bed[0]

        output_filename1 = input_filename1 + '_all_sp.txt'

        print >> sys.stderr, "Parsing alignment from %s ..." % \
                                                        (input_filename1)
        for n, exons in enumerate(parser(input_filename1), start=1):
            if len(exons) > 1:
                if parser == parse_psl:
                    add_intron(exons, db1, args.all)
                else:
                    add_intron_bed(exons, db1, args.all)

            if n % 1000 == 0:
                print >> sys.stderr, '...', n

        print >> sys.stderr, "Total introns in %s = %d\n" % \
                                            (input_filename1, len(db1))
        op = open(output_filename1, 'w')
        for junc in db1: print >> op, junc
        op.close()

    else:
        db2= set([])

        if not args.bed and not args.psl:
            raise ValueError('No input file found.')
        elif not args.bed and args.psl:
            first_file, second_file = args.psl
            first, second = adjust_parser(first_file, 'psl',
                                            second_file, 'psl',
                                            db1,
                                            db2)
        elif not args.psl and args.bed:
            first_file, second_file = args.bed
            first, second = adjust_parser(first_file, 'bed',
                                            second_file, 'bed',
                                            db1,
                                            db2)
        elif len(args.bed) == 1 and len(args.psl) == 1:
            first_file = args.bed[0]
            second_file = args.psl[0]
            first, second = adjust_parser(first_file,
                                                'bed',
                                                second_file, 'psl',
                                                db1,
                                                db2
                                            )
        elif len(args.bed) + len(args.psl) < 2:
            raise ValueError('Need two input files to compare.')

        elif len(args.bed) + len(args.psl) > 2:
            raise ValueError('Too many input files.')

        for filename, parser, db in ((first), (second)):
            print >> sys.stderr, "Parsing alignment from %s ..." % (filename)
            for n, exons in enumerate(parser(filename), start=1):
                if len(exons) > 1:
                    if parser == parse_psl:
                        add_intron(exons, db)
                    else:
                        add_intron_bed(exons, db)

                if n % 1000 == 0:
                    print >> sys.stderr, '...', n

        print >> sys.stderr, \
                    "Total introns in %s = %d\n" % (first_file, len(db1))
        print >> sys.stderr, \
                    "Total introns in %s = %d\n" % (second_file, len(db2))

        first_diff = db1.difference(db2)
        second_diff = db2.difference(db1)

        print >> sys.stderr, "Total introns not in %s = %d" % \
                                    (second_file, len(first_diff))
        print >> sys.stderr, "Total introns not in %s = %d" % \
                                    (first_file, len(second_diff))

        input_filename1 = os.path.basename(first_file)
        input_filename2 = os.path.basename(second_file)


        output_filename1 = input_filename1 + '_diff_sp.txt'
        output_filename2 = input_filename2 + '_diff_sp.txt'

        if first_diff:
            op = open(output_filename1, 'w')
            for junc in first_diff:
                print >> op, junc
            op.close()

        if second_diff:
            op = open(output_filename2, 'w')
            for junc in second_diff:
                print >> op, junc
            op.close()


if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--all', action='store_true', 
                            help='write all splice junctions to files')
    parser.add_argument('-p', '--psl', action='append', 
                            help='input file is in PSL format')
    parser.add_argument('-b', '--bed', action='append', 
                            help='input file is in BED format')

    args = parser.parse_args()
    main(args)
