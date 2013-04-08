import sys
import unittest

import networkx as nx
from utils.find_AFE import get_exon_node, Exon, find_AFE


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


class TestFindAFEPositive(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.graph = nx.DiGraph()
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
        self.graph.add_path([str(self.ex1), str(self.ex2)])
        self.transcripts = [[str(self.ex1), str(self.ex2)]]
        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)

        self.assertEqual(len(paths), 0)

    def test_positive_two_path_one_exon(self):
        path1 = [str(self.ex1), str(self.ex2)]
        path2 = [str(self.ex5), str(self.ex2)]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 2)
        self.assertItemsEqual(num_paths, [2, 2])

    def test_positive_two_path_two_exon(self):
        self.ex6 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex6)] = self.ex6
        path1 = [str(self.ex1), str(self.ex3), str(self.ex6)]
        path2 = [str(self.ex5), str(self.ex2),
                    str(self.ex3), str(self.ex6)]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]

        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 2)
        self.assertItemsEqual(num_paths, [2, 3])

    def test_positive_two_path_multiple_exon(self):
        self.ex6 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.ex7 = Exon('chrX', 9000, 10000, 'ex1.1', '+')
        self.exonsDB[str(self.ex6)] = self.ex6
        self.exonsDB[str(self.ex7)] = self.ex7

        path1 = [str(self.ex1), str(self.ex3),
                    str(self.ex6), str(self.ex7)]
        path2 = [str(self.ex5), str(self.ex7)]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]

        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 2)
        self.assertItemsEqual(num_paths, [4, 2])

    def test_positive_three_path_one_exon(self):
        path1 = [str(self.ex1), str(self.ex2)]
        path2 = [str(self.ex5), str(self.ex2)]
        path3 = [str(self.ex4), str(self.ex2)]
        self.transcripts = [path1, path2, path3]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.graph.add_path(path3)

        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 3)
        self.assertItemsEqual(num_paths, [2, 2, 2])

    def test_positive_three_path_two_exon(self):
        self.ex4 = Exon('chrX', 700, 800, 'ex1.1', '+')
        self.ex6 = Exon('chrX', 500, 600, 'ex1.1', '+')
        self.ex8 = Exon('chrX', 3500, 4500, 'ex1.1', '+')
        self.exonsDB[str(self.ex4)] = self.ex4
        self.exonsDB[str(self.ex6)] = self.ex6
        self.exonsDB[str(self.ex8)] = self.ex8

        path1 = [str(self.ex1), str(self.ex2), str(self.ex3)]
        path2 = [str(self.ex5), str(self.ex8), str(self.ex3)]
        path3 = [str(self.ex6), str(self.ex4), str(self.ex3)]
        self.transcripts = [path1, path2, path3]

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.graph.add_path(path3)

        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 3)
        self.assertItemsEqual(num_paths, [3, 3, 3])


class TestFindAFENegative(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.graph = nx.DiGraph()
        self.ex1 = Exon('chrX', 1000, 2000, 'ex1.1', '-')
        self.ex2 = Exon('chrX', 3000, 4000, 'ex1.1', '-')
        self.ex3 = Exon('chrX', 5000, 6000, 'ex1.1', '-')
        self.ex4 = Exon('chrX', 800, 900, 'ex1.1', '-')
        self.ex5 = Exon('chrX', 1500, 2500, 'ex1.1', '-')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4
        self.exonsDB[str(self.ex5)] = self.ex5

    def test_negative_one_path_one_exon(self):
        self.graph.add_path([str(self.ex1), str(self.ex2)])
        self.transcripts = [[str(self.ex1), str(self.ex2)]]
        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)

        self.assertEqual(len(paths), 0)

    def test_negative_two_path_one_exon(self):
        path1 = [str(self.ex1), str(self.ex2)]
        path2 = [str(self.ex1), str(self.ex3)]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 2)
        self.assertItemsEqual(num_paths, [2, 2])

    def test_negative_two_path_multiple_exon(self):
        self.ex6 = Exon('chrX', 7000, 8000, 'ex1.1', '-')
        self.ex7 = Exon('chrX', 9000, 10000, 'ex1.1', '-')
        self.exonsDB[str(self.ex6)] = self.ex6
        self.exonsDB[str(self.ex7)] = self.ex7

        path1 = [str(self.ex1), str(self.ex2),
                    str(self.ex3), str(self.ex6)]
        path2 = [str(self.ex1), str(self.ex7)]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]

        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 2)
        self.assertItemsEqual(num_paths, [4, 2])

    def test_negative_three_path_one_exon(self):
        path1 = [str(self.ex1), str(self.ex2)]
        path2 = [str(self.ex1), str(self.ex3)]
        path3 = [str(self.ex1), str(self.ex5)]
        self.transcripts = [path1, path2, path3]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.graph.add_path(path3)

        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 3)
        self.assertItemsEqual(num_paths, [2, 2, 2])

    def test_negative_three_path_two_exon(self):
        self.ex4 = Exon('chrX', 7500, 8500, 'ex1.1', '+')
        self.ex6 = Exon('chrX', 3500, 4500, 'ex1.1', '+')
        self.ex8 = Exon('chrX', 5500, 6500, 'ex1.1', '+')
        self.exonsDB[str(self.ex4)] = self.ex4
        self.exonsDB[str(self.ex6)] = self.ex6
        self.exonsDB[str(self.ex8)] = self.ex8

        path1 = [str(self.ex1), str(self.ex2), str(self.ex3)]
        path2 = [str(self.ex1), str(self.ex5), str(self.ex4)]
        path3 = [str(self.ex1), str(self.ex6), str(self.ex8)]
        self.transcripts = [path1, path2, path3]

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.graph.add_path(path3)

        paths = find_AFE(self.graph, self.exonsDB, self.transcripts)
        num_paths = [len(path) for path in paths]

        self.assertEqual(len(paths), 3)
        self.assertItemsEqual(num_paths, [3, 3, 3])
