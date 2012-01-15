#! /usr/local/bin/python

'''
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

parser.add_option('-o', '--output', dest='output', type='string',
                    help='output format 1=fastq, 2=fasta [default=%default]',
                    default=1, metavar='<int>')

parser.add_option('-n', '--no-duplicate', dest='noDuplicate', default=False, action='store_true', 
                    help='reads mapped to multiple locations will be \
                    written to the output only once [default=%default].')

(opts, args) = parser.parse_args()


def printFastq(name, seq, qual, outputFile):
    '''
        Print reads in FASTQ format to standard output.
    '''
    print >> outputFile, '@%s\n%s\n+\n%s' % (name, seq, qual)


def printFasta(name, seq, qual, outputFile):
    '''
        Print reads in FASTA format to standard output.
    '''
    print >> outputFile, '>%s\n%s' % (name, seq)


def writeReads(sequence, opts):
    uniques = set([])
    samfile = pysam.Samfile(args[0], 'rb')
    outputFilename = '.fq' if printOutput == printFastq else '.fa'
    with open(outputFilename, 'w') as outputFile:
        for alignedRead in samfile.fetch(sequence, opts.begin, opts.end):
            if not opts.noDuplicate:
                printOutput(alignedRead.qname,
                                alignedRead.seq,
                                alignedRead.qual,
                                outputFile)
            else:
                if alignedRead.qname not in uniques:
                    printOutput(alignedRead.qname,
                                alignedRead.seq,
                                alignedRead.qual,
                                outputFile)

            uniques.add(alignedRead.qname)


if __name__=='__main__':

    import os
    printOutput = printFastq if opts.output==1 else printFasta
    if os.path.isfile(opts.sequence):
        for line in open(opts.sequence):
            sequence = line.strip()
            writeReads(opts, sequence)
    else:
        sequence = None if opts.sequence=='all' else opts.sequence
        writeReads(opts, sequence)
