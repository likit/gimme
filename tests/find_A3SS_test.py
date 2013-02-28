import unittest
import networkx as nx

from bx.intervals import IntervalTree
from utils.find_A3SS import get_exon_node, find_A3SS, get_introns

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
    def test_find_A3SS(self):
        self.redundant = set()
        self.no_paths = []
        self.exonsDB = {}
        self.intronsDB = {}
        self.graph = nx.DiGraph()
        self.current_id = None
        self.intron_interval = IntervalTree()
        for exons, transcript_id in get_exon_node(test_file):
            self.new_id = transcript_id.split('.')[0]
            if not self.current_id: # first gene
                for e in exons: self.exonsDB[str(e)] = e
                self.graph.add_path([str(e) for e in exons])
                self.introns = get_introns(exons, self.intronsDB)

                self.assertEqual(len(exons), len(self.introns) + 1)

                for intron in self.introns:
                    self.intron_interval.insert_interval(intron)

                self.current_id = self.new_id
            else:
                if self.new_id != self.current_id:
                    if len(self.graph.nodes()) > 1:
                        self.no_paths.append(len(list(find_A3SS(self.graph,
                                                        self.intron_interval,
                                                        self.exonsDB))))
                    self.graph = nx.DiGraph()
                    self.exonsDB = {}
                    self.intronsDB = {}
                    self.current_id = self.new_id
                    self.intron_interval = IntervalTree()

                for e in exons: self.exonsDB[str(e)] = e
                self.graph.add_path([str(e) for e in exons])
                self.introns = get_introns(exons, self.intronsDB)

                self.assertEqual(len(exons), len(self.introns) + 1)

                for intron in self.introns:
                    self.intron_interval.insert_interval(intron)

        if len(self.graph.nodes()) > 1:
            self.no_paths.append(len(list(find_A3SS(self.graph,
                                                self.intron_interval,
                                                self.exonsDB))))
        self.assertEqual(self.no_paths, [0, 1, 2, 3, 4])
