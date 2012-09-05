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

        self.gene = nx.DiGraph()
        self.A = Exon('A', 'chr1', 500, 1000)
        self.B = Exon('B', 'chr1', 800, 1000)
        self.C = Exon('C', 'chr1', 2000, 2500)
        self.D = Exon('D', 'chr1', 3000, 3500)
        self.E = Exon('E', 'chr1', 4000, 4500)
        self.F = Exon('F', 'chr1', 5000, 5500)
        self.G = Exon('G', 'chr1', 5000, 5800)
        self.H = Exon('H', 'chr1', 2800, 3500)
        self.I = Exon('I', 'chr1', 3800, 4000)

        self.gene.add_edge(self.A, self.C, {'start':1000, 'end':2000})
        self.gene.add_edge(self.B, self.C, {'start':1000, 'end':2000})
        self.gene.add_edge(self.C, self.H, {'start':2500, 'end':2800})
        self.gene.add_edge(self.H, self.E, {'start':3500, 'end':4000})
        self.gene.add_edge(self.C, self.D, {'start':2500, 'end':3000})
        self.gene.add_edge(self.D, self.E, {'start':3500, 'end':4000})
        self.gene.add_edge(self.C, self.E, {'start':2500, 'end':4000})
        self.gene.add_edge(self.E, self.F, {'start':4500, 'end':5000})
        self.gene.add_edge(self.E, self.G, {'start':4500, 'end':5000})
        self.gene.add_edge(self.I, self.E, {'start':3800, 'end':4000})
        self.gene.add_edge(self.C, self.I, {'start':2500, 'end':3600})
    
    def test_graph_construction(self):
        self.assertEqual(self.gene.number_of_nodes(), 9)

    def test_identify_skipped_exon(self):
        self.assertEqual(set([(self.C, self.E)]),
                alt_iden.skipped_exon(self.gene))

    def test_identify_mutually_exclusive_exon(self):
        self.assertEqual(set([self.D,self.I]),
                alt_iden.mutual_exon(self.gene))
