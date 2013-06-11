import unittest
import networkx as nx

from utils.find_A3SS import get_exon_node, find_A3SS, Exon

'''test data contain genes with 0, 1, 2, 3 and 4 alternative splice site.

'''
test_file = "../test_data/A3SS.test.bed"


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


class TestFindA3SS(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.ex1 = Exon('chrX', 1000, 2000, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 3000, 4000, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 5000, 6000, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 3500, 4000, 'ex1.1', '+')
        self.ex5 = Exon('chrX', 1500, 2000, 'ex1.1', '+')
        self.ex6 = Exon('chrX', 3800, 4000, 'ex1.1', '+')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4
        self.exonsDB[str(self.ex5)] = self.ex5
        self.exonsDB[str(self.ex6)] = self.ex6

    def test_no_ASS(self):
        self.graph = nx.DiGraph()
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3)])
        self.graph.add_path([str(self.ex1), str(self.ex5), str(self.ex3)])
        self.events = find_A3SS(self.graph, self.exonsDB)

        self.assertEqual(len(self.events), 0)

    def test_one_ASS(self):
        self.graph = nx.DiGraph()
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3)])
        self.graph.add_path([str(self.ex1), str(self.ex4), str(self.ex3)])
        self.events = find_A3SS(self.graph, self.exonsDB)

        self.assertEqual(len(self.events), 1)
        self.assertEqual(len(self.events[0]), 2)  # two isoforms

    def test_two_ASS(self):
        self.graph = nx.DiGraph()
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3)])
        self.graph.add_path([str(self.ex1), str(self.ex4), str(self.ex3)])
        self.graph.add_path([str(self.ex1), str(self.ex6), str(self.ex3)])
        self.events = find_A3SS(self.graph, self.exonsDB)

        self.assertEqual(len(self.events), 1)
        self.assertEqual(len(self.events[0]), 3)  # three isoforms

    def test_two_events(self):
        self.ex5 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.ex6 = Exon('chrX', 9000, 10000, 'ex1.1', '+')
        self.ex7 = Exon('chrX', 7500, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex5)] = self.ex5
        self.exonsDB[str(self.ex6)] = self.ex6
        self.exonsDB[str(self.ex7)] = self.ex7
        self.graph = nx.DiGraph()
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3)])
        self.graph.add_path([str(self.ex1), str(self.ex4), str(self.ex3)])
        self.graph.add_path([str(self.ex3), str(self.ex5), str(self.ex6)])
        self.graph.add_path([str(self.ex3), str(self.ex7), str(self.ex6)])
        self.events = find_A3SS(self.graph, self.exonsDB)

        self.assertEqual(len(self.events), 2)
        self.assertEqual(len(self.events[0]), 2)  # two isoforms
        self.assertEqual(len(self.events[1]), 2)  # two isoforms
