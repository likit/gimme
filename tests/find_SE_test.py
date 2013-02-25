import unittest
import networkx as nx
from utils.find_SE import get_exon_node, find_SE

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
    def test_find_SE(self):
        self.no_paths = []
        self.graph = nx.DiGraph()
        self.current_id = None
        for exons, transcript_id in get_exon_node(test_file):
            self.new_id = transcript_id.split('.')[0]
            if not self.current_id: # first gene
                self.graph.add_path([str(e) for e in exons])
                self.current_id = self.new_id
            else:
                if self.current_id != self.new_id:
                    self.paths = find_SE(self.graph)
                    self.paths = list(find_SE(self.graph))
                    self.no_paths.append(len(self.paths))
                    self.current_id = self.new_id
                    self.graph = nx.DiGraph()

                self.graph.add_path([str(e) for e in exons])

        self.paths = list(find_SE(self.graph))
        self.no_paths.append(len(self.paths))
        self.assertEqual(self.no_paths, [1, 0, 1, 2, 3])
