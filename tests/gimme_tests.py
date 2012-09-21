'''Please run nosetests from a program main directory.'''

import sys
import os

from unittest import TestCase
import unittest
import networkx as nx

source_path = os.path.abspath('src')
if source_path not in sys.path:
    sys.path.append(os.path.abspath('src'))

import gimme

class TestCollapseExons(TestCase):
    def setUp(self):
        self.exon_db = {}
        start = 1000
        n = 1
        exons = []

        while n < 7:
            e = gimme.ExonObj('chr1', start, start + 100)
            self.exon_db[str(e)] = e
            exons.append(str(e))
            start += 300
            n += 1

        self.exon_db[exons[0]].terminal = 1 # mark a left terminal
        self.exon_db[exons[-1]].terminal = 2 # mark a right terminal

        self.exon_graph = nx.DiGraph()
        self.exon_graph.add_path(exons)

    def test_building_base_exon_db_and_exon_graph(self):
        self.assertEqual(len(self.exon_db), 6)
        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)
        self.assertItemsEqual(self.exon_graph.nodes(), ['chr1:1000-1100',
                                                    'chr1:1300-1400',
                                                    'chr1:1600-1700',
                                                    'chr1:1900-2000',
                                                    'chr1:2200-2300',
                                                    'chr1:2500-2600'])

        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600')])

    def test_collapse_left_terminal_exon(self):
        '''
            before
            L|=====|--------|======|-----R
               L|==|--------|======|-----R

            after
            L|=====|--------|======|-----R

        '''

        e = gimme.ExonObj('chr1', 1050, 1100)
        e.terminal = 1
        self.exon_db[str(e)] = e
        self.exon_graph.add_edge(str(e), 'chr1:1300-1400')
        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

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

        e = gimme.ExonObj('chr1', 2500, 2550)
        e.terminal = 2
        self.exon_db[str(e)] = e
        self.exon_graph.add_edge('chr1:2200-2300', str(e))
        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

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

        e1 = gimme.ExonObj('chr1', 700, 800)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 900, 1100)
        e2.terminal = 2
        self.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)


        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:700-800', 'chr1:900-1100'),
                    ('chr1:900-1100', 'chr1:1300-1400'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600')])

    def test_collapse_right_terminal_exon_with_other_exons(self):
        '''
            before
            L----|=====|--------|======|R
                                   L|==|-----------|===|R

            after
            L----|=====|--------|======|-----------|===|R

        '''

        e1 = gimme.ExonObj('chr1', 2550, 2600)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 2800, 2900)
        e2.terminal = 2
        self.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)


        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600'),
                    ('chr1:2500-2600', 'chr1:2800-2900'),
                    ])

    def test_collapse_right_terminal_exon_with_skipped_exon(self):
        '''
            before
            L|=====|------|=====|--------|======|R
            L|=====|---------------------|===|R

            after
            L|=====|------|=====|--------|======|R
            L|=====|---------------------|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1900, 2000)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 2500, 2550)
        e2.terminal = 2
        self.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        #print >> sys.stderr, self.exon_graph.edges()

        #nx.draw(self.exon_graph)
        #pyplot.show()
        #print >> sys.stderr, self.exon_graph.neighbors('chr1:1900-2000')

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:1900-2000', 'chr1:2500-2600'),
                    ('chr1:2200-2300', 'chr1:2500-2600'),
                    ])

    def test_collapse_left_terminal_exon_with_skipped_exon(self):
        '''
            before
              L|===|---------------------|======|R
            L|=====|------|=====|--------|======|R

            after
            L|=====|---------------------|======|R
            L|=====|------|=====|--------|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        e2.terminal = 2
        self.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1000-1100', 'chr1:1600-1700'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600'),
                    ])

    def test_collapse_left_right_terminal_exon_with_skipped_exon(self):
        '''
            before
              L|===|---------------------|======|----|===|R
            L|=====|------|=====|--------|======|----|======|R

            after
            L|=====|---------------------|======|----|======|R
            L|=====|------|=====|--------|======|----|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        self.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exons(self.exon_graph, self.exon_db)
        
        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1000-1100', 'chr1:1600-1700'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600'),
                    ])

    def test_collapse_left_utr_longer_than_100(self):
        '''
            before
                      L|========|--------|======|----|===|R
            L|=====|------|=====|--------|======|----|======|R

            after
            L|=====|------|=====|--------|======|----|======|R
                      L|========|--------|======|----|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1190, 1400)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        self.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)
        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1190-1400', 'chr1:1600-1700'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600'),
                    ])

    def test_collapse_left_utr_shorter_than_100(self):
        '''
            before
                      L|========|--------|======|----|===|R
            L|=====|------|=====|--------|======|----|======|R

            after
            L|=====|------|=====|--------|======|----|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1250, 1400)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        self.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

    def test_collapse_right_utr_longer_than_100(self):
        '''
            before
              L|===|------|=====|--------|=========|R
            L|=====|------|=====|--------|======|-------|======|R

            after
            L|=====|------|=====|--------|=========|R
            L|=====|------|=====|--------|======|-------|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1300, 1400)
        self.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1600, 1850)
        e3.terminal = 2
        self.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1300-1400', 'chr1:1600-1850'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600'),
                    ])

    def test_collapse_left_and_right_utr_longer_than_100(self):
        '''
            before
                    L|==========|--------|=========|R
            L|=====|------|=====|--------|======|-------|======|R

            after
                    L|==========|--------|=========|R
            L|=====|------|=====|--------|======|-------|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1150, 1400)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1850)
        e2.terminal = 2
        self.exon_db[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)
        self.assertItemsEqual(self.exon_graph.edges(),
                [('chr1:1000-1100', 'chr1:1300-1400'),
                    ('chr1:1300-1400', 'chr1:1600-1700'),
                    ('chr1:1150-1400', 'chr1:1600-1850'),
                    ('chr1:1600-1700', 'chr1:1900-2000'),
                    ('chr1:1900-2000', 'chr1:2200-2300'),
                    ('chr1:2200-2300', 'chr1:2500-2600'),
                    ])

    def test_collapse_left_and_right_utr_shorter_than_100(self):
        '''
            before
                    L|==========|--------|=========|R
            L|=====|------|=====|--------|======|-------|======|R

            after
            L|=====|------|=====|--------|======|-------|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1250, 1400)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1750)
        e2.terminal = 2
        self.exon_db[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

    def test_collapse_right_utr_shorter_than_100(self):
        '''
            before
              L|===|------|=====|--------|=========|R
            L|=====|------|=====|--------|======|-------|======|R

            after
            L|=====|------|=====|--------|======|-------|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1300, 1400)
        self.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1600, 1750)
        e3.terminal = 2
        self.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

    def test_collapse_all(self):
        '''
            before
                            L|==|--------|====|R
            L|=====|------|=====|--------|======|-------|======|R

            after
            L|=====|------|=====|--------|======|-------|======|R

        '''

        e1 = gimme.ExonObj('chr1', 1350, 1400)
        e1.terminal = 1
        self.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1650)
        e2.terminal = 2
        self.exon_db[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exons(self.exon_graph, self.exon_db)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

class TestAddIntrons(TestCase):
    def setUp(self):
        self.exon_db = {}
        start = 1000
        n = 1
        self.exons = []

        while n < 7:
            e = gimme.ExonObj('chr1', start, start + 100)
            self.exon_db[str(e)] = e
            self.exons.append(str(e))
            start += 300
            n += 1

        self.exon_db[self.exons[0]].terminal = 1 # mark a left terminal
        self.exon_db[self.exons[-1]].terminal = 2 # mark a right terminal

        self.intronDb = {}
        self.clusters = {}
        self.clusterNo = 0

    def test_simple(self):
        gimme.add_introns(self.exons, self.intronDb, self.exon_db,
                    self.clusters, self.clusterNo)
        self.assertEqual(len(self.intronDb), 5)


class TestMergeExons(TestCase):
    def setUp(self):
        self.exons = {'chr1':[]}
        self.e1 = gimme.ExonObj('chr1', 1000, 2000)
        self.e2 = gimme.ExonObj('chr1', 3000, 4000)
        self.e3 = gimme.ExonObj('chr1', 5000, 6000)
        self.e4 = gimme.ExonObj('chr1', 7000, 8000)


        self.exons['chr1'].append(self.e1)
        self.exons['chr1'].append(self.e2)
        self.exons['chr1'].append(self.e3)
        self.exons['chr1'].append(self.e4)

    def print_items(self, items):
        for e in items['chr1']:
            print >> sys.stderr, e,
        print >> sys.stderr, ''

    def test_no_merge_single_exons(self):
        self.exons = {'chr1':[]}
        self.exons['chr1'].append(self.e1)
        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 1)

    def test_no_merge_multiple_exons(self):
        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_subset_merge(self):
        self.e5 = gimme.ExonObj('chr1', 1100, 1800)
        self.exons['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_front(self):
        self.e5 = gimme.ExonObj('chr1', 500, 1800)
        self.exons['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_back_first(self):
        self.e5 = gimme.ExonObj('chr1', 1100, 2200)
        self.exons['chr1'].append(self.e5)
        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_back_last(self):
        self.e5 = gimme.ExonObj('chr1', 7100, 8200)
        self.exons['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exons(self.exons)
        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_back_first_last(self):
        self.e5 = gimme.ExonObj('chr1', 7100, 8200)
        self.e6 = gimme.ExonObj('chr1', 1100, 2200)
        self.exons['chr1'].append(self.e5)
        self.exons['chr1'].append(self.e6)

        self.merged_exons = gimme.merge_exons(self.exons)
        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_single_merge(self):
        self.e5 = gimme.ExonObj('chr1', 500, 8200)
        self.exons['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 1)

    def test_merge_two_exons(self):
        self.e5 = gimme.ExonObj('chr1', 1300, 3200)
        self.exons['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 3)

    def test_merge_two_exons_extend(self):
        self.e5 = gimme.ExonObj('chr1', 1300, 4200)
        self.exons['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exons(self.exons)

        self.assertEqual(len(self.merged_exons['chr1']), 3)

if __name__=='__main__':
    unittest.main()
