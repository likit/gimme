import sys
from Bio.Blast import NCBIXML


def getHits(inputFile, evalue=0.01):
    '''This function parses Blast output -- XML format
    and reports a score for each match that has e-value less
    than a given one or a default evalue (0.01).
    The max score of each query is reported to stdout.

    '''
    queryWithMatch = 0
    handle = open(inputFile)
    blastRecords = NCBIXML.parse(handle)
    for n, blastRecord in enumerate(blastRecords, start=1):
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
            queryWithMatch += 1

        if n % 100 == 0:
            print >> sys.stderr, '...', n
    pct = float(queryWithMatch) / n * 100
    print >> sys.stderr, 'Total query', n
    print >> sys.stderr, 'Total query with match %d (%.02f)' % \
                            (queryWithMatch, pct)

if __name__ == '__main__':
    inputFile = sys.argv[1]
    try:
        evalue = float(sys.argv[2])
    except IndexError:
        getHits(inputFile)
    else:
        getHits(inputFile, evalue)
