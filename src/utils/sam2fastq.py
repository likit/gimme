#! /usr/local/bin/python

'''The script extracts read alignments in BAM format and
writes them in FASTQ format. It works with both single- and
paired-end reads.

Author: Likit Preeyanon
Email: preeyano@msu.edu

'''

import pysam
from optparse import OptionParser

usage = "usage: %prog [options] bamfile"
parser = OptionParser(usage)

parser.add_option('-s', '--sequence', dest='sequence', 
                    help='a text file containing a list of reference \
                    sequences or a name of a reference sequence\
                    [default=%default]'
                    , default='all')

parser.add_option('-b', '--begin', dest='begin', type='int',
                    help='the begining of the reference sequence',
                    default=None, metavar='<int>')

parser.add_option('-e', '--end', dest='end', type='int',
                    help='the end of the reference sequence',
                    default=None, metavar='<int>')

parser.add_option('-t', '--type', dest='read_type', type='string',
                    help='a type of reads; single- or paired-end [default=%default]',
                    default='single')

parser.add_option('-f', '--output_format', dest='output_format', type='string',
                    help='output format 1=fastq, 2=fasta [default=%default]',
                    default=1, metavar='<int>')

parser.add_option('-n', '--no-duplicate', dest='noDuplicate', default=False, action='store_true', 
                    help='reads mapped to multiple locations will be \
                    written to the output only once [default=%default].')

(opts, args) = parser.parse_args()


def printFastq(alignedRead, outputFile):
    '''
        Print reads in FASTQ format to a given output file.
    '''
    print >> outputFile, '@%s\n%s\n+\n%s' % (alignedRead.qname,
                                                alignedRead.seq,
                                                alignedRead.qual)


def printFasta(alignedRead, outputFile):
    '''
        Print reads in FASTA format to a given output file.
    '''
    print >> outputFile, '>%s\n%s' % (alignedRead.qname, alignedRead.seq)

def printFastqPaired(alignedRead, mate, outputFile1, outputFile2):
    '''
        Print reads in FASTQ format to a given output file.
    '''
    print >> outputFile1, '@%s\n%s\n+\n%s' % (alignedRead.qname,
                                                alignedRead.seq,
                                                alignedRead.qual)

    print >> outputFile2, '@%s\n%s\n+\n%s' % (mate.qname,
                                                mate.seq,
                                                mate.qual)


def printFastaPaired(alignedRead, mate, outputFile1, outputFile2):
    '''
        Print reads in FASTA format to a given output file.
    '''
    print >> outputFile1, '>%s\n%s' % (alignedRead.qname, alignedRead.seq)
    print >> outputFile2, '>%s\n%s' % (mate.qname, mate.seq)

def writeReads(sequence, opts):
    uniques = set([])

    samfile = pysam.Samfile(args[0], 'rb')
    if opts.read_type == 'single':
        outputFilename = sequence+'.fq' if printOutput == printFastq else sequence+'.fa'

        with open(outputFilename, 'w') as outputFile:
            for n, alignedRead in enumerate(samfile.fetch(sequence, opts.begin, opts.end), start=1):
                if not opts.noDuplicate:
                    printOutput(alignedRead, outputFile)
                else:
                    if alignedRead.qname not in uniques:
                        printOutput(alignedRead.qname, outputFile)

                uniques.add(alignedRead.qname)
                if n % 1000 == 0:
                    print >> sys.stderr, '...', sequence, n

    elif opts.read_type == 'paired':
        outputFilename1 = sequence+'_1.fq' if printOutput == printFastqPaired else sequence+'_1.fa'
        outputFilename2 = sequence+'_2.fq' if printOutput == printFastqPaired else sequence+'_2.fa'

        outputFile1 = open(outputFilename1, 'w')
        outputFile2 = open(outputFilename2, 'w')
        for n, alignedRead in enumerate(samfile.fetch(sequence, opts.begin, opts.end), start=1):
            if not alignedRead.is_unmapped and alignedRead.is_read1:
                if alignedRead.is_proper_pair:
                    if not alignedRead.mate_is_unmapped and \
                            alignedRead.mate_is_reverse:
                        if not opts.noDuplicate:
                            mate = samfile.mate(alignedRead)
                            printOutput(alignedRead,
                                            mate,
                                            outputFile1,
                                            outputFile2)
                        else:
                            if alignedRead.qname not in uniques:
                                mate = samfile.mate(alignedRead)
                                printOutput(alignedRead,
                                            mate,
                                            outputFile1,
                                            outputFile2)

                        uniques.add(alignedRead.qname)
            if n % 1000 == 0:
                print >> sys.stderr, '...', sequence, n
        outputFile1.close()
        outputFile2.close()

    else:
        raise ValueError, 'Unsupported read type.'


if __name__=='__main__':

    import os, sys
    
    if opts.read_type == 'single':
        printOutput = printFastq if opts.output_format==1 else printFasta
    elif opts.read_type == 'paired':
        printOutput = printFastqPaired if opts.output_format==1 else printFastaPaired
    else:
        raise ValueError, 'Unsupported read type.'

    if os.path.isfile(opts.sequence):
        for line in open(opts.sequence):
            sequence = line.strip()
            writeReads(sequence, opts)
            print >> sys.stderr, '...' + sequence + ' done.'
    else:
        sequence = None if opts.sequence=='all' else opts.sequence
        writeReads(sequence, opts)
