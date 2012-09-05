'''This script identifies several kinds of alternative splicings
from given gene models.

The script now only supports gene models in BED format.

'''

import networkx as nx

def skipped_exon(gene):
    '''Identify skipped exons, if exists in a gene.
    gene is a DiGraph object containing exons and introns as
    nodes and edges respectively.

    '''

    skipped_exons = set([])
    for node in gene.nodes():
        neighbors = gene.neighbors(node)
        for nb in neighbors:
            for s in gene.successors(nb):
                if s in neighbors:
                    skipped_exons.add((node,s))

    return skipped_exons


def mutual_exon(gene):
    '''Identify mutually exclusive exons in a given gene.'''
    for node in gene.nodes():
        if gene.out_degree(node) > 1:
            candidates = {}
            for s in gene.successors(node):
                candidates[s.name] = set([])
                for ss in nx.algorithms.dfs_successors(gene, s):
                    if ss != s:
                        candidates[s.name].add(ss.name)
                print node.name, candidates
                clusters = {}
            for exon,de in candidates.iteritems():
                clusters[exon] = []
                for ex,d in candidates.iteritems():
                    if d.intersection(de):
                        clusters[exon].append(ex)
            print clusters
