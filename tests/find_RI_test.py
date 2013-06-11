import unittest
import networkx as nx

from utils.find_RI import get_exon_node, find_RI, add_intervals, Exon
from bx.intervals import IntervalTree

'''test data contain genes with 0, 1, 2, 3 and 4 retained introns.'''
test_file = "../test_data/RI.test.bed"


class TestLoadData(unittest.TestCase):
    def test_load_data(self):
        for exons, transcript_id in get_exon_node(test_file):
            pass
        self.assertTrue(len(exons) > 1)


class TestExonGraph(unittest.TestCase):
    def test_build_exon_graph(self):
        for exons, transcript_id in get_exon_node(test_file):
            self.graph = nx.DiGraph()
            self.graph.add_path([str(e) for e in exons])
            self.assertEqual(len(exons), len(self.graph.nodes()))


class TestAddIntervals(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.ex1 = Exon('chrX', 1000, 2000, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 3000, 4000, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 5000, 6000, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4
        self.graph = nx.DiGraph()
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex3), str(self.ex4)])

        self.tree = add_intervals(self.graph, self.exonsDB)

        self.assertEqual(len(self.tree.find(1500, 5500)), 3)


class TestFindRI(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.ex1 = Exon('chrX', 1000, 2000, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 3000, 4000, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 5000, 6000, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4
        self.tree = IntervalTree()
        self.tree.add_interval(self.ex1)
        self.tree.add_interval(self.ex2)
        self.tree.add_interval(self.ex3)
        self.tree.add_interval(self.ex4)
        self.graph = nx.DiGraph()

    def test_no_retained_introns(self):
        self.path1 = [str(self.ex1), str(self.ex2), str(self.ex3)]
        self.path2 = [str(self.ex1), str(self.ex3), str(self.ex4)]
        self.graph.add_path(self.path1)
        self.graph.add_path(self.path2)
        self.events = list(find_RI(self.graph, self.tree, self.exonsDB))

        self.assertEqual(len(self.events), 0)

    def test_one_retained_introns(self):
        self.ex5 = Exon('chrX', 3000, 6000, 'ex1.1', '+')
        self.exonsDB[str(self.ex5)] = self.ex5
        self.tree.add_interval(self.ex5)

        self.path1 = [str(self.ex1), str(self.ex2),
                        str(self.ex3), str(self.ex4)]
        self.path2 = [str(self.ex1), str(self.ex5), str(self.ex4)]
        self.graph.add_path(self.path1)
        self.graph.add_path(self.path2)
        self.events = list(find_RI(self.graph, self.tree, self.exonsDB))

        self.assertEqual(len(self.events), 1)

    def test_two_retained_introns(self):
        self.ex5 = Exon('chrX', 1000, 4000, 'ex1.1', '+')
        self.exonsDB[str(self.ex5)] = self.ex5
        self.tree.add_interval(self.ex5)

        self.ex6 = Exon('chrX', 5000, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex6)] = self.ex6
        self.tree.add_interval(self.ex6)

        self.path1 = [str(self.ex1), str(self.ex2),
                        str(self.ex3), str(self.ex4)]
        self.path2 = [str(self.ex5), str(self.ex6)]
        self.graph.add_path(self.path1)
        self.graph.add_path(self.path2)
        self.events = list(find_RI(self.graph, self.tree, self.exonsDB))

        self.assertEqual(len(self.events), 2)
