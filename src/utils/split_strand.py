import sys
import string

import networkx as nx
from pygr import seqdb

genome = seqdb.SequenceFileDB('../chick_4.fa')
table = string.maketrans('ACGT', 'TGCA')

def get_splice_sites(exon1, exon2):
    chrom, pos = exon1.split(':')
    start, end = [int(p) for p in pos.split('-')]
    donor = genome[chrom][end:end + 2]

    chrom, pos = exon2.split(':')
    start, end = [int(p) for p in pos.split('-')]
    acceptor = genome[chrom][start - 2:start]

    return str(donor), str(acceptor)

def identify_strand(splice_sites):
    donor, acceptor = splice_sites

    canonicals = set([('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')])
    rev_canonicals = set([('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')])
    if (donor, acceptor) in canonicals:
        return 1
    elif (donor, acceptor) in rev_canonicals:
            return -1
    else:
        return 0

def compare_edges(edge):
    return int(edge[0].split(':')[1].split('-')[0])

def split(graph):

    class Edgeobj(object):
        def __init__(self, edge, ss, strand):
            self.edge = edge
            self.ss = ss
            self.strand = strand

    # print 'total nodes %d' % graph.number_of_nodes()
    edges = {}

    pos_graph = nx.DiGraph(strand='+')  # a graph for positive strand
    neg_graph = nx.DiGraph(strand='-')  # a graph for negative strand

    strand_scores = []
    sorted_edges = sorted(graph.edges(), key=compare_edges)
    for edge in sorted_edges:
        splice_sites = get_splice_sites(*edge)
        strand = identify_strand(splice_sites)
        edges[edge] = Edgeobj(edge, splice_sites, strand)
        strand_scores.append(strand)

    score_matrix = [sum(strand_scores[0:3]) / 3.0]
    i = 1
    while (i < len(sorted_edges) - 1):
        score_matrix.append(sum(strand_scores[i-1:i+2]) / 3.0)
        i += 1
    score_matrix.append(sum(strand_scores[-3:]) / 3.0)

    for i in range(len(sorted_edges)):
        edge = edges[sorted_edges[i]]
        print >> sys.stderr, \
                sorted_edges[i], edge.ss, strand_scores[i], score_matrix[i]

    for i in range(len(sorted_edges)):
        if score_matrix[i] > 0:
            pos_graph.add_edge(*sorted_edges[i])
        elif score_matrix[i] < 0:
            neg_graph.add_edge(*sorted_edges[i])

    neg_graph.remove_edges_from(pos_graph.edges())
    pos_graph.remove_edges_from(neg_graph.edges())

    return pos_graph, neg_graph
