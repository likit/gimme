'''The script reads gene models in BED format and searches for
minimum isoforms.

If the total number in original gene models is lesser than that
of minimum isoforms, the script return the original transcripts.

'''

import sys, csv

import networkx as nx

class ExonObj(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)

def parseBed(filename):
    '''Reads BED file and returns exons of a transcript.'''

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

            yield geneId, exons, row

def create_bipartite_graph(G):
    '''Return a bipartite graph with top nodes = G.nodes and
    bottom nodes = {1,2,3...} which 1=G.nodes[1], etc.
    Each edge has capacity=1.0.

    G is a directed graph. Start and End nodes are not included
    in a bipartite graph.

    '''
    B = nx.DiGraph()
    g = G.copy()
    g.remove_nodes_from(['Start', 'End'])

    bottom_nodes = {}
    node_index = [] # note, nodes are not id'd in topological order

    node_id = 0
    for node in g.nodes():
        B.add_edge('S', node, capacity=1.0) # add edge to [S]ource node

        bottom_nodes[node] = node_id
        node_index.append(node)

        B.add_edge(node_id, 'T', capacity=1.0) # add edge to [T]arget node

        node_id += 1
    # print bottom_nodes;

    for node in g.nodes():
        for e in g[node].keys():
            node_id = bottom_nodes[e]
            B.add_edge(node, node_id, capacity=1.0)

    # for edge in B.edges():
    #     print edge
    return B, node_index

def run_max_flow(G):
    '''Return edges in max flow analysis.

    G is a bipartite graph created from create_bipartite_graph().

    '''
    edges = nx.max_flow_min_cost(G, 'S', 'T')
    return edges

def rebuild_edges(edges, node_index):
    '''Rebuild edges from bipartite graph to orginal directed graph.

    edges are output from run_max_flow().

    '''
    K = nx.DiGraph()
    for node in edges:
        if node == 'S':
            continue # ignore S node
        for e in edges[node]:
            if e == 'T':
                continue # ignore T node

            # print node, e

            if edges[node][e] > 0.0:
                K.add_edge(node, node_index[e])

    return K.edges()

def remove_maxflow_edges(mf_edges, B):
    '''Remove edges from max flow path of bipartite graph B.'''

    for node in mf_edges:
        if node == 'S': continue # ignore S node
        for e in mf_edges[node]:
            if e == 'T':
                continue # ignore T node
            if mf_edges[node][e] == 1.0:
                B[node][e]['capacity'] = 0.0 # disable the edge

def add_path(edges, paths, G):
    '''Build paths from a given edges.

    G is an original directed graph.

    '''
    K = nx.DiGraph()
    J = nx.DiGraph()
    K.add_edges_from(edges)
    for SG in nx.algorithms.weakly_connected_component_subgraphs(K):
        head_node = None
        tail_node = None
        for node in SG.nodes():
            if not SG.predecessors(node):
                head_node = node
            if not SG.successors(node):
                tail_node = node

        # print 'head node= %s, tail node = %s' % (head_node, tail_node)

        body = nx.algorithms.shortest_path(SG, head_node, tail_node)
        head = nx.algorithms.shortest_path(G, 'Start', head_node)
        tail = nx.algorithms.shortest_path(G, tail_node, 'End')

        SG.add_path(head)
        SG.add_path(body)
        SG.add_path(tail)

        path = [p for p in nx.algorithms.dfs_preorder_nodes(SG, 'Start')]

        J.add_path(path)
        if set(J.edges()).difference(set(G.edges())):
            raise ValueError, "Error: edges are added."
        # for e in SG.edges():
        #     print e
        path_str = '->'.join(path)
        paths.add(path_str)

def get_min_paths(G, verbose=True):
    '''Returns minimal paths including all edges.
    G is a directed graph.

    '''
    total_edges = len(G.edges())
    remaining_edges = total_edges
    paths = set() # store unique paths
    B, node_index = create_bipartite_graph(G)

    mf_edges = run_max_flow(B)
    edges = set(rebuild_edges(mf_edges, node_index))

    mf_round = 1
    while edges: # a max flow path is found
        remaining_edges -= len(edges)
        if (total_edges > 50) and verbose: # display progress
            if verbose:
                print >> sys.stderr, \
                    '\t... #%d found %d edges' % (mf_round, len(edges))

        add_path(edges, paths, G)
        remove_maxflow_edges(mf_edges, B)
        mf_edges = run_max_flow(B)
        edges = set(rebuild_edges(mf_edges, node_index))
        mf_round += 1

    g = G.copy()
    g.remove_nodes_from(['Start', 'End'])

    paths = list(paths)
    for i in range(len(paths)):
        paths[i] = paths[i].split('->')
        paths[i].remove('Start')
        paths[i].remove('End')

    K = nx.DiGraph()
    for path in paths:
        K.add_path(path)

    if set(K.edges()) != set(g.edges()):
        raise ValueError, "Error: Some edges are added or removed."

    return paths

def printBed(path, db, geneId, transId):
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

def make_graph(bedfile, exon_db):
    def get_path(exons, exon_db):
        path = [str(e) for e in exons]
        for exon in exons:
            exon_db[str(exon)] = exon
        path.insert(0, 'Start')
        path.append('End')

        return path

    gene_id = None
    total = 0
    G = nx.DiGraph()
    transcripts = []
    for id, exons, trns in parseBed(bedfile):
        if not gene_id:
            gene_id = id

        if gene_id == id:
            path = get_path(exons, exon_db)
            G.add_path(path)
            total += 1
            transcripts.append(trns)
        else:
            yield gene_id, G, total, transcripts
            transcripts = []
            total = 1
            G = nx.DiGraph()
            gene_id = id
            path = get_path(exons, exon_db)
            G.add_path(path)

    yield gene_id, G, total, transcripts

def main(argv, verbose=True):
    bedfile = argv[1]
    exon_db = {}
    for n, items in enumerate(make_graph(bedfile, exon_db), start=1):
        gene_id, G, total, transcripts = items

        if verbose:
            print >> sys.stderr, '%s, total nodes = %d\n\tSearch' %\
                                            (gene_id, len(G.nodes()))

        paths = get_min_paths(G, verbose)
        if len(G.edges()) < total:
            for n, path in enumerate(paths, start=1):
                printBed(path, exon_db, gene_id, n)
            if verbose:
                print >> sys.stderr, \
                    '\ttotal max = %d, min = %d' % (total, len(paths))
        else:
            for trns in transcripts:
                print '\t'.join(trns)
            if verbose:
                print >> sys.stderr, \
                    '\ttotal original = %d' % total

if __name__=='__main__':
    main(sys.argv)
    # G = nx.DiGraph()
    # G.add_path(['Start', 'A', 'B', 'C', 'E', 'F', 'G', 'I', 'J', 'End'])
    # G.add_path(['Start', 'A', 'B', 'C', 'E', 'F', 'G', 'H', 'End'])
    # paths = get_min_paths(G)
    # print len(paths)
