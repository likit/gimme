import sys, os

source_path = os.path.abspath('src')
if source_path not in sys.path:
    sys.path.append(os.path.abspath('src'))

import networkx as nx
from unittest import TestCase

from gimme import ExonObj, collapseExons

class TestCollapseExons(TestCase):
    def setUp(self):
        self.exonDb = {}
        start = 100
        n = 1
        exons = []

        while n < 7:
            e = ExonObj('chr1', start, start + 100)
            self.exonDb[str(e)] = e
            exons.append(str(e))
            start += 300
            n += 1

        self.exonDb[exons[0]].terminal = 1 # mark a left terminal
        self.exonDb[exons[-1]].terminal = 2 # mark a right terminal

        self.exon_graph = nx.DiGraph()
        self.exon_graph.add_path(exons)

    def test_building_base_exon_db_and_exon_graph(self):
        self.assertEqual(len(self.exonDb), 6)
        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

    def test_collapse_left_terminal_exon(self):
        '''
            before
            L|=====|--------|======|-----R
               L|==|--------|======|-----R

            after
            L|=====|--------|======|-----R

        '''

        e = ExonObj('chr1', 120, 200)
        e.terminal = 1
        self.exonDb[str(e)] = e
        self.exon_graph.add_edge(str(e), 'chr1:400-500')
        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

    def test_collapse_right_terminal_exon(self):
        '''
            before
            L----|=====|--------|======|R
            L----|=====|--------|===|R

            after
            L----|=====|--------|======|R

        '''

        e = ExonObj('chr1', 1600, 1650)
        e.terminal = 2
        self.exonDb[str(e)] = e
        self.exon_graph.add_edge('chr1:1300-1400', str(e))
        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

    def test_collapse_left_terminal_exon_with_other_exons(self):
        '''
            before
                         |===|--------|======|R
            L|====|----|=====|R

            after
            L|====|----|=====|--------|======|R

        '''

        e1 = ExonObj('chr1', 10, 50)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 70, 200)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)


        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

    def test_collapse_right_terminal_exon_with_other_exons(self):
        '''
            before
            L----|=====|--------|======|R
                                   L|==|--------|===|R

            after
            L----|=====|--------|======|-----------|===|R

        '''

        e1 = ExonObj('chr1', 1650, 1700)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1900, 2000)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)


        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)
