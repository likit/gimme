import unittest
import networkx as nx
from utils.find_SE import get_exon_node, find_SE, Exon

'''test data contain genes with 0, 1, 2 and 3 exon skipping.'''
test_file = "../test_data/SE.test.bed"


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


class TestFindSE(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.graph = nx.DiGraph()
        self.ex1 = Exon('chrX', 1000, 2000, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 3000, 4000, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 5000, 6000, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4

    def test_no_skipped_exon(self):
        self.ex5 = Exon('chrX', 9000, 10000, 'ex1.1', '+')
        self.exonsDB[str(self.ex5)] = self.ex5
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3)])
        self.graph.add_path([str(self.ex3), str(self.ex4), str(self.ex5)])
        self.events = list(find_SE(self.graph))

        self.assertEqual(len(self.events), 0)

    def test_one_skipped_exon(self):
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3)])
        self.graph.add_path([str(self.ex1), str(self.ex3)])
        self.events = list(find_SE(self.graph))

        self.assertEqual(len(self.events), 1)  # one skipped exon
        self.assertEqual(len(self.events[0]) + 1, 2)  # two isoforms

    def test_one_double_skipped_exon(self):
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex3), str(self.ex4)])
        self.graph.add_path([str(self.ex1), str(self.ex4)])
        self.events = list(find_SE(self.graph))

        self.assertEqual(len(self.events), 1)  # one skipped exon
        self.assertEqual(len(self.events[0]) + 1, 2)  # two isoforms

    def test_two_single_skipped_exon(self):
        self.ex5 = Exon('chrX', 9000, 10000, 'ex1.1', '+')
        self.exonsDB[str(self.ex5)] = self.ex5

        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3)])
        self.graph.add_path([str(self.ex1), str(self.ex3)])

        self.graph.add_path([str(self.ex3), str(self.ex4), str(self.ex5)])
        self.graph.add_path([str(self.ex3), str(self.ex5)])
        self.events = list(find_SE(self.graph))

        self.assertEqual(len(self.events), 2)  # one skipped exon
        self.assertEqual(len(self.events[0]) + 1, 2)  # two isoforms
        self.assertEqual(len(self.events[1]) + 1, 2)  # two isoforms
