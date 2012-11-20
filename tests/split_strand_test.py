import networkx as nx
import utils.split_strand as sp
from unittest import TestCase

class TestMakeGraph(TestCase):
    def test_simple_digraph(self):
        self.G = nx.DiGraph()
        self.G.add_path([1,2,3,4,5,6,7])
        self.assertEqual(7, len(self.G.nodes()))

    def test_simple_digraph_single_strand(self):
        self.G = nx.DiGraph()
        self.G.add_path([1,2,3,4,5,6,7], strand='+')
        self.assertEqual('+', self.G[2][3]['strand'])

    def test_simple_digraph_double_strand(self):
        self.G = nx.DiGraph()
        self.G.add_path([1,2,3,4,5,6,7], strand='+')
        self.G.add_path([7,6,5,4], strand='-')
        self.assertEqual('+', self.G[2][3]['strand'])
        self.assertEqual('-', self.G[6][5]['strand'])

class TestSplitLinearGraph(TestCase):
    def setUp(self):
        self.G = nx.DiGraph()
        self.G.add_path([1,2,3,4,5,6,7], strand='+')
        self.G.add_path([4,5,6,7], strand='-')

    def test_simple_two_strands(self):
        strands, g1, g2 = sp.split(self.G)
        self.assertEqual(2, len(strands)) # number of strands

    # def test_simple_two_strands_with_unknown(self):
    #     self.G.add_edge(4,5, strand=None)
    #     strands, g1, g2 = sp.split(self.G)
    #     self.assertEqual(3, len(strands)) # number of strands

    def test_simple_split(self):
        strands, g1, g2 = sp.split(self.G)
        self.assertEqual(4, g1.number_of_nodes())
        self.assertEqual(4, g2.number_of_nodes())

    # def test_split_one_mismatch(self):
    #     self.G.add_path([4,5,6,7], strand='-')
    #     self.G.add_edge(5,6, strand='+')
    #     strands, g1, g2 = sp.split(self.G)
    #     self.assertEqual(4, g1.number_of_nodes())
    #     self.assertEqual(4, g2.number_of_nodes())
