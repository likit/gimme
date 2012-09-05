import sys


def parseBed(filename):
    for line in open(filename):
        cols = line.split()
        name = cols[3]
        strand = cols[5]
        yield name, strand


def assignStrand(junctionFile, transcripts):
    for line in open(junctionFile):
        name, _, ts = line.strip().split()
        tlist = ts.split(',')
        strands = []
        for t in tlist:
            strands.append(transcripts[t])

        if strands.count('+') > strands.count('-'):
            print '%s\t%s' % (name, '+')
        else:
            print '%s\t%s' % (name, '-')


def main(bedFile, junctionFile):
    transcripts = {}
    for name, strand in parseBed(bedFile):
        transcripts[name] = strand
    
    assignStrand(junctionFile, transcripts)


if __name__=='__main__':
    bedFile = sys.argv[1]
    junctionFile = sys.argv[2]
    main(bedFile, junctionFile)
