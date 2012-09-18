''' The program searches for a minimum set of paths in a graph
that contains all edges. It is designed to work with gene models
in BED format. The output is written to stdout.

Author: Likit Preeyanon
Email : preeyano@msu.edu

'''

import sys, csv
from sys import stderr

import networkx as nx


class ExonObj(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


def parseBed(filename):
    '''Reads BED file and returns a group of exons
    of a transcript at a time.

    '''

    with open(filename) as fp:
        for row in csv.reader(fp, dialect='excel-tab'):
            exons = []
            chrom = row[0]
            chromStart = int(row[1])
            geneId = row[3].split('.')[0]

            '''Get all exons except terminal ones.'''
            blockStarts = [int(i) for i in row[-1].split(',')]
            blockSizes = [int(i) for i in row[-2].split(',')]

            if len(blockStarts) == 1:
                continue

            for i in range(len(blockStarts)):
                start = chromStart + blockStarts[i]
                end = start + blockSizes[i]
                exons.append(ExonObj(chrom, start, end))

            yield geneId, exons


def getEdges(exons):
    n = nx.DiGraph()
    n.add_path([str(e) for e in exons])
    return n.edges()


def walkEdges(edge, paths, pathCounts,
                    subPaths, subPathEdges, uniqEdges, ignored, l):

    '''A recursive function that walks through different
    edges between two given paths. Edges get added to sub-path
    edges. The fuction returns when all edges get added.

    '''

    ranked = sorted([e for e in uniqEdges[edge]],
                                key=lambda x: pathCounts[x], reverse=True)

    if len(ranked) > 1:
        '''Remove all paths that have been visited.'''

        newRanked = [e for e in ranked if e not in ignored]

        '''If there is no edge left, do nothing.
        Otherwise, asign a new list to ranked.

        '''
        if newRanked:
            ranked = newRanked

    #print >> stderr, edge, ranked, 'ignored', ignored,

    subPaths.add(ranked[0])
    subPathEdges.update(set(paths[ranked[0]]))
    #print >> stderr, subPaths, len(subPathEdges), len(uniqEdges)

    if subPathEdges == set(uniqEdges):
        return

    for p in ranked:
        diffEdges = set(paths[p]).difference(set(paths[ranked[0]]))

        #print >> stderr, '-'*l, 'move to', p

        if diffEdges:
            ignored.add(p)
            for edge in diffEdges:
                if subPaths.intersection(uniqEdges[edge]):
                    continue

                walkEdges(edge, paths, pathCounts,
                        subPaths, subPathEdges, uniqEdges, ignored, l + 2)


def getMinPaths(paths):
    uniqEdges = {}
    for i in range(len(paths)):
        for edge in paths[i]:
            try:
                uniqEdges[edge].append(i)
            except KeyError:
                uniqEdges[edge] = [i]

    pathCounts = {}
    for e in uniqEdges.itervalues():
        for path in e:
            try:
                pathCounts[path] += 1
            except KeyError:
                pathCounts[path] = 1
    '''
    #-------- for debugging purposes-------.
    for
    sortedUniqEdges = sorted(uniqEdges, key=lambda x: len(uniqEdges[x]))
    for e in sortedUniqEdges:
        ranked = sorted(uniqEdges[e], key=lambda x: x)
        print >> stderr, e, ranked,
    '''
    

    subPathEdges = set([])
    subPaths = set([])
    ignored = set([])
    sortedUniqEdges = sorted(uniqEdges, key=lambda x: len(uniqEdges[x]))
    for edge in sortedUniqEdges:
        walkEdges(edge, paths, pathCounts,
                        subPaths, subPathEdges, uniqEdges, ignored, 1)
        if subPathEdges == set(uniqEdges):
            break

    
    assert subPathEdges == set(uniqEdges), \
                            'Cannot get all unique edges. %d %d' %\
                            (len(subPathEdges), len(uniqEdges))
    if __name__=='__main__':
        print >> stderr, ', %d minimal isoforms' % (len(subPaths))

    minPaths = []
    for p in subPaths:
        g = nx.DiGraph()
        g.add_edges_from(paths[p])
        minPaths.append(g.nodes())

    return minPaths


def printBED(path, db, geneId, transId):
    path = sorted(path, key=lambda x: db[x].start)
    firstExon = db[path[0]]
    lastExon = db[path[-1]]
    chrom = firstExon.chrom
    chromStart = firstExon.start
    chromEnd = lastExon.end

    blockStarts = [str(db[e].start - chromStart) for e in path]
    blockSizes = [str(db[e].end - db[e].start) for e in path]
    blockCount = len(blockStarts)
    
    name = '%s.%d' % (geneId, transId)
    strand = '+'
    score = 1000
    thickStart = chromStart
    thickEnd = chromEnd
    itemRgb='0,0,0'

    writer = csv.writer(sys.stdout, dialect='excel-tab')

    writer.writerow((chrom,
                    chromStart,
                    chromEnd,
                    name,
                    score,
                    strand,
                    thickStart,
                    thickEnd,
                    itemRgb,
                    blockCount,
                    ','.join(blockSizes),
                    ','.join(blockStarts)))
def main(filename):
    db = {}
    currGeneId = None
    print >> sys.stderr, 'processing gene models from %s ...' % filename
    for n, (geneId, exons) in enumerate(parseBed(filename)):
        db.update(dict([(str(exon), exon) for exon in exons]))

        if currGeneId == None:
            currGeneId = geneId
            paths = []
            paths.append(getEdges(exons))
        elif currGeneId == geneId:
            paths.append(getEdges(exons))
        else:
            transId = 0
            print >> sys.stderr, '%s - %d maximal isoforms' %\
                                        (currGeneId, len(paths)),
            minPaths = getMinPaths(paths)
            for path in minPaths:
                printBED(path, db, currGeneId, transId)
                transId += 1

            paths = [getEdges(exons)]
            currGeneId = geneId

    # Process the last transcript
    transId = 0
    print >> sys.stderr, '%s - %d maximal isoforms' %\
                                (currGeneId, len(paths)),
    minPaths = getMinPaths(paths)
    for path in minPaths:
        printBED(path, db, currGeneId, transId)
        transId += 1

if __name__=='__main__':
    inputFile = sys.argv[1]
    main(inputFile)
