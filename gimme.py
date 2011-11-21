import sys
import csv

import networkx as nx
import psl_parser
import matplotlib.pyplot as plot

from sys import stderr, stdout

MIN_INTRON = 21
MIN_UTR = 100

class ExonObj:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.terminal = None
        self.nextExons = set([])
        self.clusters = []
        self.introns = set([])

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


def parsePSL(filename):
    '''Reads alignments from PSL format and create
    exon objects from each transcript.

    '''

    for pslObj in psl_parser.read(open(filename)):
        exons = []
        for i in range(len(pslObj.attrib['tStarts'])):
            exonStart = pslObj.attrib['tStarts'][i]
            exonEnd = exonStart + pslObj.attrib['blockSizes'][i]

            exon = ExonObj(pslObj.attrib['tName'], exonStart, exonEnd)
            exons.append(exon)

        yield exons


def addIntrons(exons, intronDb, exonDb, clusters, clusterNo):
    '''Get introns from a set of exons.'''

    existingClusters = set([])

    introns = []

    for i in range(len(exons)):
        currExon = exonDb[str(exons[i])]
        try:
            nextExon = exonDb[str(exons[i + 1])]
        except IndexError:
            pass
        else:
            currExon.nextExons.add(str(nextExon))

            intronStart = currExon.end + 1
            intronEnd = nextExon.start - 1

            intron = ExonObj(currExon.chrom, intronStart, intronEnd)

            try:
                intron_ = intronDb[str(intron)]
            except KeyError:

                intronDb[str(intron)] = intron
                intron.exons = set([str(currExon), str(nextExon)])
                introns.append(intron)
            else:
                introns.append(intron_)
                intron_.exons.add(str(currExon))
                intron_.exons.add(str(nextExon))
                existingClusters.add(intron_.cluster)

            currExon.introns.add(str(intron))
            nextExon.introns.add(str(intron))

    assert len(introns) == len(exons) - 1

    if not existingClusters:
        g = nx.DiGraph()
        if len(introns) > 1:
            g.add_path([str(i) for i in introns])
            for intron in introns:
                intron.cluster = clusterNo
                for exon in intron.exons:
                    exonDb[exon].clusters.append(clusterNo)
        else:
            g.add_node(str(introns[0]))
            introns[0].cluster = clusterNo
            for exon in introns[0].exons:
                exonDb[exon].clusters.append(clusterNo)

    else:
        g = nx.DiGraph(exons=set([]))
        for cl in existingClusters:
            g.add_edges_from(clusters[cl].edges())
            k = clusters.pop(cl)

        for intron in g.nodes():
            intronDb[intron].cluster = clusterNo
            for exon in intronDb[intron].exons:
                exonDb[exon].clusters.append(clusterNo)

        if len(introns) > 1:
            g.add_path([str(i) for i in introns])

            for intron in introns:
                intron.cluster = clusterNo
                for exon in intron.exons:
                    exonDb[exon].clusters.append(clusterNo)
        else:
            g.add_node(str(introns[0]))
            introns[0].cluster = clusterNo
            for exon in introns[0].exons:
                exonDb[exon].clusters.append(clusterNo)

    clusters[clusterNo] = g

    return clusterNo


def collapseExons(g, exonDb, intronDb):
    mergedExons = set([])
    exons = [exonDb[e] for e in g.nodes()]
    sortedExons = sorted(exons, key=lambda x: (x.end, x.start))

    i = 0
    currExon = sortedExons[i]
    while i <= len(sortedExons):
        try:
            nextExon = sortedExons[i + 1]
        except IndexError:
            pass
        else:
            if currExon.end == nextExon.end:
                if nextExon.terminal == 1:
                    g.add_edges_from([(str(currExon), n)\
                            for n in g.neighbors(str(nextExon))])
                    g.remove_node(str(nextExon))
                    mergedExons.add(str(nextExon))
                else:
                    if currExon.terminal == 1:
                        if nextExon.start - currExon.start <= MIN_UTR:
                            g.add_edges_from([(str(nextExon), n)\
                                    for n in g.neighbors(str(currExon))])
                            g.remove_node(str(currExon))
                            mergedExons.add(str(currExon))

                    currExon = nextExon
            else:
                currExon = nextExon
        i += 1

    i = 0
    exons = [exonDb[e] for e in g.nodes()]
    sortedExons = sorted(exons, key=lambda x: (x.start, x.end))
    currExon = sortedExons[0]
    while i <= len(sortedExons):
        try:
            nextExon = sortedExons[i + 1]
        except IndexError:
            pass
        else:
            if currExon.start == nextExon.start:
                if currExon.terminal == 2:
                    g.add_edges_from([(str(nextExon), n)\
                            for n in g.neighbors(str(currExon))])

                    g.remove_node(str(currExon))
                    mergedExons.add(str(currExon))
                    currExon = nextExon
                else:
                    if nextExon.terminal == 2:
                        if nextExon.end - currExon.end <= MIN_UTR:
                            g.add_edges_from([(str(currExon), n)\
                                    for n in g.neighbors(str(nextExon))])

                            g.remove_node(str(nextExon))
                            mergedExons.add(str(nextExon))
                        else:
                            currExon = nextExon
                    else:
                        currExon = nextExon
            else:
                currExon = nextExon
        i += 1
    return mergedExons


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
        i += 1

    newExon.append(currExon)
    return newExon


def addExon(db, exons):
    '''Adds exons to a given dictionary.
    db      : a dictionary to which exons will be added.
    exons   : a list of exon objects.

    '''

    '''Change a terminal attribute of a left end exon
    and a right end exon to True. For a middile exon,
    a terminal attribute has value 'None' by default.
    
    '''
    exons[0].terminal = 1
    exons[-1].terminal = 2
    if str(exons[0]) == "chrZ:23487400-23487466":
        print >> stderr, [str(exon) for exon in exons]
        print >> stderr, '%s:%d-%d' % (exons[0].chrom, exons[0].end + 1, exons[-1].start - 1)

    for exon in exons:
        try:
            exon_ = db[str(exon)]
        except KeyError:
            db[str(exon)] = exon
        else:
            if not exon.terminal and exon_.terminal:
                exon_.terminal = None


def walkDown(intronCoord, path, allPath, cluster):
    '''Returns all downstream exons from a given exon.
    
    '''

    if cluster.successors(intronCoord) == []:
        allPath.append(path[:])
        return
    else:
        for nex in cluster.successors(intronCoord):
            if nex not in path:
                path.append(nex)

            walkDown(nex, path, allPath, cluster)

            path.pop()


def getPath(cluster):
    '''Returns all paths of a given cluster.'''

    roots = [node for node in cluster.nodes() \
                    if not cluster.predecessors(node)]
    allPaths = []

    for root in roots:
        path = [root]
        walkDown(root, path, allPaths, cluster)

    return allPaths


def buildSpliceGraph(cluster, intronDb, exonDb, mergedExons):
    allPaths = getPath(cluster)

    g = nx.DiGraph()
    for path in allPaths:
        for intronCoord in path:
            intron = intronDb[intronCoord]
            for exon in intron.exons:
                for nextExon in exonDb[exon].nextExons:
                    if nextExon not in mergedExons:
                        g.add_edge(exon, nextExon)

    return g


def printBedGraph(transcript, geneId, tranId):
    '''Print a splice graph in BED format.'''
    exons = sorted([exonDb[e] for e in transcript],
                            key=lambda x: (x.start, x.end))

    chromStart = exons[0].start
    chromEnd = exons[-1].end
    chrom = exons[0].chrom

    blockStarts = ','.join([str(exon.start - chromStart) for exon in exons])
    blockSizes = ','.join([str(exon.end - exon.start) for exon in exons])

    name = '%s:%d.%d' % (chrom, geneId, tranId)
    score = 1000
    itemRgb = '0,0,0'
    thickStart = chromStart
    thickEnd = chromEnd
    strand = '+'
    blockCount = len(exons)

    writer = csv.writer(stdout, dialect='excel-tab')
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
                    blockSizes,
                    blockStarts))


def buildGeneModels(exonDb, intronDb, clusters):
    print >> stderr, 'Building gene models...'
    removedClusters = set([])
    passedExons = set([])
    numTranscripts = 0
    geneId = 0

    for n, exon in enumerate(exonDb.itervalues(), start=1):
        if exon.clusters and str(exon) not in passedExons:
            g = nx.DiGraph()
            for cl in set(exon.clusters):
                if cl not in removedClusters and cl in clusters:
                    for intron in clusters[cl].nodes():
                        for exon in intronDb[intron].exons:
                            for nextExon in exonDb[exon].nextExons:
                                if nextExon not in passedExons:
                                    g.add_edge(exon, nextExon)

                    removedClusters.add(cl)

            if g.nodes():
                geneId += 1
                transId = 1
                passedExons = passedExons.union(collapseExons(g, exonDb, intronDb))

                for transcript in getPath(g):
                    printBedGraph(transcript, geneId, transId)
                    numTranscripts += 1
                    transId += 1

                for exon in g.nodes():
                    passedExons.add(exon)

        if n % 1000 == 0:
            print >> stderr, '...', n
    return geneId, numTranscripts


if __name__=='__main__':

    inputFile = sys.argv[1]

    exonDb = {}
    intronDb = {}
    clusters = {}
    clusterNo = 0

    print >> stderr, 'Parsing alignments from %s...' % inputFile
    for n, exons in enumerate(parsePSL(inputFile), start=1):
        exons = deleteGap(exons)  # delete small intron (<= MIN_INTRON)

        if len(exons) > 1:
            addExon(exonDb, exons)
            addIntrons(exons, intronDb, exonDb, clusters, clusterNo)
            clusterNo += 1
        else:
            exons[0].single = True

        if n % 1000 == 0:
            print >> stderr, '...', n

    geneId, numTranscripts = buildGeneModels(exonDb, intronDb, clusters)
    print >> stderr, '\nTotal exons = %d' % len(exonDb)
    print >> stderr, 'Total genes = %d' % geneId
    print >> stderr, 'Total transcripts = %d' % (numTranscripts)
    print >> stderr, 'Isoform/gene = %.2f' % (float(numTranscripts) / len(clusters))
