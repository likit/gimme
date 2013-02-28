import unittest
import networkx as nx

from utils.find_RI import get_exon_node, find_RI, add_intervals

'''test data contain genes with 0, 1, 2 and 3 exon skipping.'''
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


class TestFindSE(unittest.TestCase):
    def test_find_RI(self):
        self.exonsDB = {}
        self.graph = nx.DiGraph()
        self.current_id = None
        self.no_paths = []
        for exons, transcript_id in get_exon_node(test_file):
            self.new_id = transcript_id.split('.')[0]
            # print >> sys.stderr, current_id, new_id
            if not self.current_id:  # first gene
                for e in exons:
                    self.exonsDB[str(e)] = e
                self.graph.add_path([str(e) for e in exons])
                self.current_id = self.new_id
            else:
                if self.new_id != self.current_id:
                    self.interval_tree = add_intervals(self.graph,
                                                        self.exonsDB)
                    self.no_paths.append(len(list(find_RI(self.graph,
                                                        self.interval_tree,
                                                        self.exonsDB))))

                    self.graph = nx.DiGraph()
                    self.exonsDB = {}
                    self.current_id = self.new_id

                for e in exons:
                    self.exonsDB[str(e)] = e
                self.graph.add_path([str(e) for e in exons])

        self.interval_tree = add_intervals(self.graph, self.exonsDB)
        self.no_paths.append(len(list(find_RI(self.graph,
                                            self.interval_tree,
                                            self.exonsDB))))

        self.assertEqual(self.no_paths, [1, 2, 3, 4])
