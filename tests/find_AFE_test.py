import sys
import unittest

import networkx as nx
from utils.find_AFE import get_exon_node, add_path, Exon, find_AFE


# class TestLoadData(unittest.TestCase):
#     def test_load_data(self):
#         for exons, transcript_id in get_exon_node(test_file):
#             pass
#         self.assertTrue(len(exons) > 1)


# class TestExonGraph(unittest.TestCase):
#     def test_build_exon_graph(self):
#         for exons, transcript_id in get_exon_node(test_file):
#             self.graph = nx.DiGraph()
#             self.graph.add_path([str(e) for e in exons])
#             self.assertEqual(len(exons), len(self.graph.nodes()))


class TestAddPath(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.graph = nx.Graph()
        self.ex1 = Exon('chrX', 1000, 2000, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 3000, 4000, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 5000, 6000, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 7000, 8000, 'ex1.1', '+')

    def test_positive_strand(self):
        exons = [self.ex1, self.ex2, self.ex3, self.ex4]
        add_path(exons, self.graph)
        self.assertItemsEqual(self.graph.nodes(),
                                            ['start',
                                                'chrX:1000-2000',
                                                'chrX:3000-4000',
                                                'end',
                                            ])

    def test_negative_strand(self):
        self.ex1.strand = '-'
        self.ex2.strand = '-'
        self.ex3.strand = '-'
        self.ex4.strand = '-'

        exons = [self.ex1, self.ex2, self.ex3, self.ex4]
        add_path(exons, self.graph)
        self.assertItemsEqual(self.graph.nodes(),
                                            ['start',
                                                'chrX:7000-8000',
                                                'chrX:5000-6000',
                                                'end',
                                            ])


class TestFindAFE(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.graph = nx.Graph()
        self.ex1 = Exon('chrX', 1000, 2000, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 3000, 4000, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 5000, 6000, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 800, 900, 'ex1.1', '+')
        self.ex5 = Exon('chrX', 1500, 2500, 'ex1.1', '+')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4
        self.exonsDB[str(self.ex5)] = self.ex5

    def test_positive_one_path_one_exon(self):
        self.graph.add_path(['start', str(self.ex1),
                                str(self.ex2), 'end'])
        paths = list(find_AFE(self.graph))
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 1)
        self.assertItemsEqual(num_paths, [2])

    def test_positive_two_path_one_exon(self):
        self.graph.add_path(['start', str(self.ex1),
                                str(self.ex2), 'end'])
        self.graph.add_path(['start', str(self.ex5),
                                str(self.ex2), 'end'])
        paths = list(find_AFE(self.graph))
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 2)
        self.assertItemsEqual(num_paths, [2, 2])

    def test_positive_two_path_two_exon(self):
        self.graph.add_path(['start', str(self.ex1),
                                str(self.ex3), 'end'])
        self.graph.add_path(['start', str(self.ex5),
                                str(self.ex2), str(self.ex3), 'end'])
        paths = list(find_AFE(self.graph))
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 2)
        self.assertItemsEqual(num_paths, [2, 3])

    def test_positive_three_path_one_exon(self):
        self.graph.add_path(['start', str(self.ex1),
                                str(self.ex2), 'end'])
        self.graph.add_path(['start', str(self.ex5),
                                str(self.ex2), 'end'])
        self.graph.add_path(['start', str(self.ex4),
                                str(self.ex2), 'end'])

        paths = list(find_AFE(self.graph))
        # num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 3)

    def test_positive_three_path_two_exon(self):
        self.ex4 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.graph.add_path(['start', str(self.ex1),
                                str(self.ex2), str(self.ex4), 'end'])
        self.graph.add_path(['start', str(self.ex5),
                                str(self.ex2), str(self.ex4), 'end'])
        self.graph.add_path(['start', str(self.ex2),
                                str(self.ex4), 'end'])

        paths = list(find_AFE(self.graph))
        # num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 3)
