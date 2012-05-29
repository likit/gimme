''' This script reports differential splice junctions 
of two gene models. It can read models in PSL and BED format.
Please adjust the parser function according to the data format.

All splice junctions from input files are written in first_sp.txt and
second_sp.txt respectively.

Differential introns are written in first_diff_sp.txt and
second_diff_sp.txt.
These files will not be created if there is no differential splice
junctions.

Author: Likit Preeyanon
Email: preeyano@msu.edu

'''


import csv, os.path, sys
from collections import namedtuple

import pslparser

MIN_INTRON = 21

Transcript = namedtuple('Transcript', ['chrom',
                                        'start',
                                        'end',
                                        'blockSizes',
                                        'blockStarts',
                                        'name',
                                        ])

class Exon(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.setName()
        self.introns = set([])
        self.terminal = None
        self.nextExons = set([])
        self.prevExons = set([])
    
    def setName(self):
        self.name = self.__str__()

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


def deleteGap(exons):
    i = 0

    newExon = []
    currExon = exons[i]

    while True:
        try:
            nextExon = exons[i + 1]
        except IndexError:
            break
        else:
            if nextExon.start - currExon.end <= MIN_INTRON:
                currExon.end = nextExon.end
            else:
                newExon.append(currExon)
                currExon = nextExon

            currExon.setName()
        i += 1

    newExon.append(currExon)
    return newExon


def parsePSL(filename):
    '''Reads alignments from PSL format and create
    exon objects from each transcript.

    '''

    for pslObj in pslparser.read(open(filename)):
        exons = []
        for i in range(len(pslObj.attrib['tStarts'])):
            exonStart = pslObj.attrib['tStarts'][i]
            exonEnd = exonStart + pslObj.attrib['blockSizes'][i]

            exon = Exon(pslObj.attrib['tName'], exonStart, exonEnd)
            exons.append(exon)

        yield exons


def addIntron(exons, intronDb):
    for i in range(len(exons)):
        currExon = exons[i]
        try:
            nextExon = exons[i + 1]
        except IndexError:
            continue
        else:
            intron = "%s:%d-%d" % (currExon.chrom,
                                        currExon.end,
                                        nextExon.start)

            try:
                assert currExon.end < nextExon.start, '%s %s %s' % (currExon, nextExon, intron)
            except AssertionError:
                print >> sys.stderr, '%s skipped ' % intron
            else:
                if nextExon.start - currExon.end > MIN_INTRON:
                    intronDb.add(intron)


def parseBED(bedFile):

    with open(bedFile) as fp:
        for row in csv.reader(fp, dialect='excel-tab'):
            blockStarts = [int(i) for i in row[11].split(',')]
            blockSizes = [int(i) for i in row[10].split(',')]

            transcript = Transcript(row[0],
                                    int(row[1]),
                                    int(row[2]),
                                    blockSizes,
                                    blockStarts,
                                    row[3],
                                    )
            yield transcript


def addIntronBed(transcript, modelIntronDb):
    exons = []

    for i in range(len(transcript.blockStarts)):
        start = transcript.blockStarts[i] + transcript.start
        end = transcript.blockSizes[i] + start
        exon = Exon(transcript.chrom, start, end)

        exons.append(exon)

    addIntron(exons, modelIntronDb)


def adjustParser(inputFile1, format1, inputFile2, format2, db1, db2):
    parsers = {'psl':parsePSL, 'bed':parseBED}

    parser1 = parsers[format1]
    parser2 = parsers[format2]
    first = (inputFile1, parser1, db1)
    second = (inputFile2, parser2, db2)

    return first, second

def main(args):
    db1 = set([])
    if args.all:
        '''Only print out all junctions without comparing junctions.

        Only one input file needed.

        '''
        if args.psl:
            parser = parsePSL
            inputFileName1 = args.psl[0]
        elif args.bed:
            parser = parseBED
            inputFileName1 = args.bed[0]
        else: # default format
            parser = parseBED
            inputFileName1 = args.bed[0]

        outputFileName1 = inputFileName1 + '_all_sp.txt'

        print >> sys.stderr, "Parsing alignment from %s ..." % (inputFileName1)
        for n, exons in enumerate(parser(inputFileName1), start=1):

            if len(exons) > 1:
                '''Skip a single exon.'''

                if parser == parsePSL:
                    addIntron(exons, db1)
                else:
                    addIntronBed(exons, db1)

            if n % 1000 == 0:
                print >> sys.stderr, '...', n

        print >> sys.stderr, "Total introns in %s = %d\n" % (inputFileName1, len(db1))

        op = open(outputFileName1, 'w')
        for junc in db1:
            print >> op, junc
        op.close()

    else:
        db2= set([])

        if not args.bed and not args.psl:
            raise ValueError, 'No input file found.'
        elif not args.bed and args.psl:
            firstFile, secondFile = args.psl
            first, second = adjustParser(firstFile, 'psl',
                                            secondFile, 'psl',
                                            db1,
                                            db2)
        elif not args.psl and args.bed:
            firstFile, secondFile = args.bed
            first, second = adjustParser(firstFile, 'bed',
                                            secondFile, 'bed',
                                            db1,
                                            db2)
        elif len(args.bed) == 1 and len(args.psl) == 1:
            firstFile = args.bed[0]
            secondFile = args.psl[0]
            first, second = adjustParser(firstFile, 'bed',
                                            secondFile, 'psl',
                                            db1,
                                            db2)
        elif len(args.bed) + len(args.psl) < 2:
            raise ValueError, 'Need two input files to compare.'

        elif len(args.bed) + len(args.psl) > 2:
            raise ValueError, 'Too many input files.'

        for filename, parser, db in ((first), (second)):
            print >> sys.stderr, "Parsing alignment from %s ..." % (filename)
            for n, exons in enumerate(parser(filename), start=1):

                if len(exons) > 1:
                    '''Skip a single exon.'''

                    if parser == parsePSL:
                        addIntron(exons, db)
                    else:
                        addIntronBed(exons, db)

                if n % 1000 == 0:
                    print >> sys.stderr, '...', n

        print >> sys.stderr, "Total introns in %s = %d\n" % (firstFile, len(db1))
        print >> sys.stderr, "Total introns in %s = %d\n" % (secondFile, len(db2))

        firstDiff = db1.difference(db2)
        secondDiff = db2.difference(db1)

        print >> sys.stderr, "Total introns not in %s = %d" % \
                                    (secondFile, len(firstDiff))
        print >> sys.stderr, "Total introns not in %s = %d" % \
                                    (firstFile, len(secondDiff))

        inputFileName1 = os.path.basename(firstFile)
        inputFileName2 = os.path.basename(secondFile)


        outputFileName1 = inputFileName1 + '_diff_sp.txt'
        outputFileName2 = inputFileName2 + '_diff_sp.txt'

        if firstDiff:
            op = open(outputFileName1, 'w')
            for junc in firstDiff:
                print >> op, junc
            op.close()

        if secondDiff:
            op = open(outputFileName2, 'w')
            for junc in secondDiff:
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
