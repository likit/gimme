'''This script identifies several kinds of alternative splicings
from given gene models.

The script now only supports gene models in BED format.

'''

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
