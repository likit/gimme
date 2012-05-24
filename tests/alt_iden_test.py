import unittest
import networkx as nx
from collections import namedtuple

from gimme.src.utils import alt_iden

Exon = namedtuple('Exon', ['name', 'chrom', 'start', 'end'])

class test_graph(unittest.TestCase):
    def setUp(self):
        '''The graph contains one skipped exon,
        an alternative splice site and two alternative UTRs.

        '''

        self.G = nx.DiGraph()
        A = Exon('A', 'chr1', 500, 1000)
        B = Exon('B', 'chr1', 800, 1000)
        C = Exon('C', 'chr1', 2000, 2500)
        D = Exon('D', 'chr1', 3000, 3500)
        E = Exon('E', 'chr1', 4000, 4500)
        F = Exon('F', 'chr1', 5000, 5500)
        G = Exon('G', 'chr1', 5000, 5800)
        H = Exon('H', 'chr1', 2800, 3500)
        I = Exon('I', 'chr1', 3800, 4000)

        self.G.add_edge(A, C, {'start':1000, 'end':2000})
        self.G.add_edge(B, C, {'start':1000, 'end':2000})
        self.G.add_edge(C, H, {'start':2500, 'end':2800})
        self.G.add_edge(H, E, {'start':3500, 'end':4000})
        self.G.add_edge(C, D, {'start':2500, 'end':3000})
        self.G.add_edge(D, E, {'start':3500, 'end':4000})
        self.G.add_edge(C, E, {'start':2500, 'end':4000})
        self.G.add_edge(E, F, {'start':4500, 'end':5000})
        self.G.add_edge(E, G, {'start':4500, 'end':5000})
        self.G.add_edge(I, E, {'start':3800, 'end':4000})
        self.G.add_edge(C, I, {'start':2500, 'end':3600})
    
    def test_graph_construction(self):
        self.assertEqual(self.G.number_of_nodes(), 9)

    def test_identify_skipped_exon(self):
        pass
        self.assertEqual(set([('C', 'E')]), alt_iden.skipped_exon(self.G))
