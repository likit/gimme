import sys
from Bio.Blast import NCBIXML

def getHits(inputFile, evalue=0.01):
    '''This function parses Blast output -- XML format
    
    and calculate bitscore/length ratio for each match
    that has e-value less than a given one. The max score
    of each queryNames is reported to stdout.

    '''
    handle = open(inputFile)
    blastRecords = NCBIXML.parse(handle)
    for blastRecord in blastRecords:
        queryName = blastRecord.query.split()[0]
        maxScore = None
        maxScoreSubj = None
        for alignment in blastRecord.alignments:
            for hsp in alignment.hsps:
                subject = alignment.title.split()[0]
                if hsp.expect < evalue:
                    if hsp.bits > maxScore:
                        maxScore = hsp.bits
                        maxScoreSubj = subject
        if maxScore:
            print >> sys.stdout, '%s\t%.16f\t%s' % (queryName,
                                                    maxScore,
                                                    maxScoreSubj,
                                                    )
            print >> sys.stderr, '%s\t%.16f\t%s' % (queryName,
                                                    maxScore,
                                                    maxScoreSubj,
                                                    )

if __name__=='__main__':
    inputFile = sys.argv[1]
    try:
        evalue = float(sys.argv[2])
    except IndexError:
        getHits(inputFile)
    else:
        getHits(inputFile, evalue)
