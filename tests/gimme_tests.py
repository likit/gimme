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
        self.align_db = gimme.AlignmentDB()
        start = 1000
        n = 1
        exons = []

        while n < 7:
            e = gimme.ExonObj('chr1', start, start + 100)
            self.align_db.exon_db[str(e)] = e
            exons.append(str(e))
            start += 300
            n += 1

        self.align_db.exon_db[exons[0]].terminal = 1 # mark a left terminal
        self.align_db.exon_db[exons[-1]].terminal = 2 # mark a right terminal

        self.exon_graph = nx.DiGraph()
        self.exon_graph.add_path(exons)

    def test_building_base_exon_db_and_exon_graph(self):
        self.assertEqual(len(self.align_db.exon_db), 6)
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
        self.align_db.exon_db[str(e)] = e
        self.exon_graph.add_edge(str(e), 'chr1:1300-1400')
        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e)] = e
        self.exon_graph.add_edge('chr1:2200-2300', str(e))
        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 900, 1100)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)


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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 2800, 2900)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)


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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 2500, 2550)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.align_db.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exon(self.exon_graph, self.align_db)
        
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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.align_db.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1700)
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.align_db.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1300, 1400)
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1600, 1850)
        e3.terminal = 2
        self.align_db.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1850)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1750)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1300, 1400)
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 1600, 1750)
        e3.terminal = 2
        self.align_db.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exon(self.exon_graph, self.align_db)

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
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1600, 1650)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        gimme.collapse_exon(self.exon_graph, self.align_db)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

    def test_two_no_collapse_one_collapse(self):
        '''
            before
          L|=========|------|=====|R
            L|=|-|=|--------|=====|R
           L|========|------|=====|--------|======|-------|======|R

            after
            L|=|-|=|--------|=====|R
          L|=========|------|=====|--------|======|-------|======|R

        '''
        e1 = gimme.ExonObj('chr1', 990, 1010)
        e1.terminal = 1
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 1300, 1400)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 950, 1100)
        e3.terminal = 1
        self.align_db.exon_db[str(e3)] = e3

        e4 = gimme.ExonObj('chr1', 1030, 1040)
        self.align_db.exon_db[str(e4)] = e4

        self.exon_graph.add_edge(str(e1), str(e4))
        self.exon_graph.add_edge(str(e3), str(e2))
        self.exon_graph.add_edge(str(e4), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 9)
        self.assertEqual(len(self.exon_graph.edges()), 8)

        gimme.collapse_exon(self.exon_graph, self.align_db)

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

    def test_no_collapse2(self):
        '''
            before
            L|=====|------|===========|R
                          L|=|-|======|------|====|R
            after
            L|=====|------|===========|R
                          L|=|-|======|------|====|R

        '''
        e1 = gimme.ExonObj('chr1', 2510, 2520)
        e1.terminal = 1
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 2530, 2600)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 2700, 2800)
        e3.terminal = 1
        self.align_db.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 9)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exon(self.exon_graph, self.align_db)

        self.assertEqual(len(self.exon_graph.nodes()), 9)
        self.assertEqual(len(self.exon_graph.edges()), 7)

    def test_no_collapse3(self):
        '''
            before
            L|=====|------|======|R
                     L|=|-|======|------|====|R
            after
            L|=====|------|======|R
                     L|=|-|======|------|====|R

        '''
        e1 = gimme.ExonObj('chr1', 2400, 2450)
        e1.terminal = 1
        self.align_db.exon_db[str(e1)] = e1

        e2 = gimme.ExonObj('chr1', 2500, 2600)
        e2.terminal = 2
        self.align_db.exon_db[str(e2)] = e2

        e3 = gimme.ExonObj('chr1', 2700, 2800)
        e3.terminal = 1
        self.align_db.exon_db[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        gimme.collapse_exon(self.exon_graph, self.align_db)

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

class TestAddIntrons(TestCase):
    def setUp(self):
        self.align_db = gimme.AlignmentDB()
        start = 1000
        n = 1
        self.exons = []

        while n < 7:
            e = gimme.ExonObj('chr1', start, start + 100)
            self.exons.append(e)
            start += 300
            n += 1

        self.exons[0].terminal = 1 # mark a left terminal
        self.exons[-1].terminal = 2 # mark a right terminal

    def test_simple(self):
        gimme.add_intron(self.exons, self.align_db, {}, 0)

        self.assertEqual(len(self.align_db.intron_db), 5)


class TestMergeExons(TestCase):
    def setUp(self):
        self.align_db = gimme.AlignmentDB()
        self.align_db.single_exons_db= {'chr1':[]}
        self.e1 = gimme.ExonObj('chr1', 1000, 2000)
        self.e2 = gimme.ExonObj('chr1', 3000, 4000)
        self.e3 = gimme.ExonObj('chr1', 5000, 6000)
        self.e4 = gimme.ExonObj('chr1', 7000, 8000)


        self.align_db.single_exons_db['chr1'].append(self.e1)
        self.align_db.single_exons_db['chr1'].append(self.e2)
        self.align_db.single_exons_db['chr1'].append(self.e3)
        self.align_db.single_exons_db['chr1'].append(self.e4)

    def print_items(self, items):
        for e in items['chr1']:
            print >> sys.stderr, e,
        print >> sys.stderr, ''

    def test_no_merge_single_exons(self):
        self.align_db.single_exons_db = {'chr1':[]}
        self.align_db.single_exons_db['chr1'].append(self.e1)
        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 1)

    def test_no_merge_multiple_exons(self):
        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_subset_merge(self):
        self.e5 = gimme.ExonObj('chr1', 1100, 1800)
        self.align_db.single_exons_db['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_front(self):
        self.e5 = gimme.ExonObj('chr1', 500, 1800)
        self.align_db.single_exons_db['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_back_first(self):
        self.e5 = gimme.ExonObj('chr1', 1100, 2200)
        self.align_db.single_exons_db['chr1'].append(self.e5)
        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_back_last(self):
        self.e5 = gimme.ExonObj('chr1', 7100, 8200)
        self.align_db.single_exons_db['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exon(self.align_db)
        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_extend_back_first_last(self):
        self.e5 = gimme.ExonObj('chr1', 7100, 8200)
        self.e6 = gimme.ExonObj('chr1', 1100, 2200)
        self.align_db.single_exons_db['chr1'].append(self.e5)
        self.align_db.single_exons_db['chr1'].append(self.e6)

        self.merged_exons = gimme.merge_exon(self.align_db)
        self.assertEqual(len(self.merged_exons['chr1']), 4)

    def test_single_merge(self):
        self.e5 = gimme.ExonObj('chr1', 500, 8200)
        self.align_db.single_exons_db['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 1)

    def test_merge_two_exons(self):
        self.e5 = gimme.ExonObj('chr1', 1300, 3200)
        self.align_db.single_exons_db['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 3)

    def test_merge_two_exons_extend(self):
        self.e5 = gimme.ExonObj('chr1', 1300, 4200)
        self.align_db.single_exons_db['chr1'].append(self.e5)

        self.merged_exons = gimme.merge_exon(self.align_db)

        self.assertEqual(len(self.merged_exons['chr1']), 3)


class TestSplitExonGroups(TestCase):
    max_intron = 200

    def setUp(self):
        self.align_db = gimme.AlignmentDB()
        self.exons = []
        e1 = gimme.ExonObj('chr1', 1000, 1100)
        e2 = gimme.ExonObj('chr1', 1200, 1300)

        self.exons.append(e1)
        self.exons.append(e2)

    def test_no_split_one_exon(self):
        self.exons = []
        e1 = gimme.ExonObj('chr1', 1000, 1100)

        self.exons.append(e1)

        split = gimme.remove_large_intron(self.exons,
                                            TestSplitExonGroups.max_intron,
                                        )
        self.assertEqual(len(split), 1)

    def test_no_split_two_exons(self):
        split = gimme.remove_large_intron(self.exons,
                                            TestSplitExonGroups.max_intron,
                                        )
        self.assertEqual(len(split), 1)

    def test_split_one_two_exons(self):
        self.exons = []
        e1 = gimme.ExonObj('chr1', 1000, 1100)
        e2 = gimme.ExonObj('chr1', 2800, 2900)

        self.exons.append(e1)
        self.exons.append(e2)

        split = gimme.remove_large_intron(self.exons,
                                            TestSplitExonGroups.max_intron,
                                        )
        self.assertEqual(len(split), 2)

    def test_split_one_back(self):
        e3 = gimme.ExonObj('chr1', 1800, 1900)

        self.exons.append(e3)

        split = gimme.remove_large_intron(self.exons,
                                            TestSplitExonGroups.max_intron,
                                        )
        self.assertEqual(len(split), 2)

    def test_split_one_front(self):
        e3 = gimme.ExonObj('chr1', 500, 700)

        self.exons.insert(0, e3)

        split = gimme.remove_large_intron(self.exons,
                                            TestSplitExonGroups.max_intron,
                                        )
        self.assertEqual(len(split), 2)

    def test_split_two(self):
        e3 = gimme.ExonObj('chr1', 1800, 1900)
        e4 = gimme.ExonObj('chr1', 2100, 2200)
        e5 = gimme.ExonObj('chr1', 3000, 3400)
        e6 = gimme.ExonObj('chr1', 3500, 3600)

        self.exons += [e3, e4, e5, e6]

        split = gimme.remove_large_intron(self.exons,
                                            TestSplitExonGroups.max_intron,
                                        )
        self.assertEqual(len(split), 3)

if __name__=='__main__':
    unittest.main()
