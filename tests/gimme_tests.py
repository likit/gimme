'''Please run nosetests from a program main directory.'''

import sys
import os
import StringIO
import csv

from unittest import TestCase
import networkx as nx

source_path = os.path.abspath('src')
if source_path not in sys.path:
    sys.path.append(os.path.abspath('src'))

from gimme import ExonObj, collapseExons, addIntrons, parsePSL

class TestCollapseExons(TestCase):
    def setUp(self):
        self.exonDb = {}
        start = 1000
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

        e = ExonObj('chr1', 1050, 1100)
        e.terminal = 1
        self.exonDb[str(e)] = e
        self.exon_graph.add_edge(str(e), 'chr1:1300-1400')
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

        e = ExonObj('chr1', 2500, 2550)
        e.terminal = 2
        self.exonDb[str(e)] = e
        self.exon_graph.add_edge('chr1:2200-2300', str(e))
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

        e1 = ExonObj('chr1', 700, 800)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 900, 1100)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)


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

        e1 = ExonObj('chr1', 2550, 2600)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 2800, 2900)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)


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

        e1 = ExonObj('chr1', 1900, 2000)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 2500, 2550)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1600, 1700)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2
        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 7)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1600, 1700)
        self.exonDb[str(e2)] = e2

        e3 = ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.exonDb[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        collapseExons(self.exon_graph, self.exonDb)
        
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

        e1 = ExonObj('chr1', 1190, 1400)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1600, 1700)
        self.exonDb[str(e2)] = e2

        e3 = ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.exonDb[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1250, 1400)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1600, 1700)
        self.exonDb[str(e2)] = e2

        e3 = ExonObj('chr1', 1900, 1950)
        e3.terminal = 2
        self.exonDb[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1300, 1400)
        self.exonDb[str(e2)] = e2

        e3 = ExonObj('chr1', 1600, 1850)
        e3.terminal = 2
        self.exonDb[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1150, 1400)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1600, 1850)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1250, 1400)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1600, 1750)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1050, 1100)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1300, 1400)
        self.exonDb[str(e2)] = e2

        e3 = ExonObj('chr1', 1600, 1750)
        e3.terminal = 2
        self.exonDb[str(e3)] = e3

        self.exon_graph.add_edge(str(e1), str(e2))
        self.exon_graph.add_edge(str(e2), str(e3))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 7)

        collapseExons(self.exon_graph, self.exonDb)

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

        e1 = ExonObj('chr1', 1350, 1400)
        e1.terminal = 1
        self.exonDb[str(e1)] = e1

        e2 = ExonObj('chr1', 1600, 1650)
        e2.terminal = 2
        self.exonDb[str(e2)] = e2

        self.exon_graph.add_edge(str(e1), str(e2))

        self.assertEqual(len(self.exon_graph.nodes()), 8)
        self.assertEqual(len(self.exon_graph.edges()), 6)

        collapseExons(self.exon_graph, self.exonDb)

        self.assertEqual(len(self.exon_graph.nodes()), 6)
        self.assertEqual(len(self.exon_graph.edges()), 5)

class TestAddIntrons(TestCase):
    def setUp(self):
        self.exonDb = {}
        start = 1000
        n = 1
        self.exons = []

        while n < 7:
            e = ExonObj('chr1', start, start + 100)
            self.exonDb[str(e)] = e
            self.exons.append(str(e))
            start += 300
            n += 1

        self.exonDb[self.exons[0]].terminal = 1 # mark a left terminal
        self.exonDb[self.exons[-1]].terminal = 2 # mark a right terminal

        self.intronDb = {}
        self.clusters = {}
        self.clusterNo = 0

    def test_simple(self):
        addIntrons(self.exons, self.intronDb, self.exonDb,
                    self.clusters, self.clusterNo)
        self.assertEqual(len(self.intronDb), 5)

class TestParseExons(TestCase):
    def setUp(self):
        self.input_file = StringIO.StringIO()

    def test_simple(self):
        fake_data = [
                0, # matches
                0, # misMatches
                0, # repMatches
                3, # nCount
                9, # qNumInsert
                9, # qBaseInsert
                19,# tNumInsert
                9, # tBaseInsert
                '+', # strand
                'transcript', # query sequence name
                300, # qSize
                0, # qStart
                300, # qEnd
                'chr1', # tName
                3000, # tSize
                1000, # tStart
                1900, # tEnd
                4, # blockCount
                '100,100,100,100,', # blockSize, the last comma is needed
                '0,300,600,800,', # qStarts, the last comma is needed
                '0,300,600,800,', # tStarts, the last comma is needed
                ]
        writer = csv.writer(self.input_file, dialect='excel-tab')
        writer.writerow(fake_data)
        self.input_file.seek(0)

        selected_exons = parsePSL(self.input_file, 20000).next()
        self.assertEqual(len(selected_exons), 1)
        self.assertEqual(len(selected_exons[0]), 4)

    def test_remove_one_large_intron_one_exon_set(self):
        fake_data = [
                0, # matches
                0, # misMatches
                0, # repMatches
                3, # nCount
                9, # qNumInsert
                9, # qBaseInsert
                19,# tNumInsert
                9, # tBaseInsert
                '+', # strand
                'transcript', # query sequence name
                300, # qSize
                0, # qStart
                300, # qEnd
                'chr1', # tName
                3000, # tSize
                1000, # tStart
                31000, # tEnd
                4, # blockCount
                '100,100,100,100,', # blockSize, the last comma is needed
                '0,300,600,30000,', # qStarts, the last comma is needed
                '0,300,600,30000,', # tStarts, the last comma is needed
                ]
        writer = csv.writer(self.input_file, dialect='excel-tab')
        writer.writerow(fake_data)
        self.input_file.seek(0)

        selected_exons = parsePSL(self.input_file, 20000).next()
        self.assertEqual(len(selected_exons), 1)
        self.assertEqual(len(selected_exons[0]), 3)

    def test_remove_one_large_intron_two_exon_sets(self):
        fake_data = [
                0, # matches
                0, # misMatches
                0, # repMatches
                3, # nCount
                9, # qNumInsert
                9, # qBaseInsert
                19,# tNumInsert
                9, # tBaseInsert
                '+', # strand
                'transcript', # query sequence name
                300, # qSize
                0, # qStart
                300, # qEnd
                'chr1', # tName
                3000, # tSize
                1000, # tStart
                30400, # tEnd
                5, # blockCount
                '100,100,100,100,100,', # blockSize, the last comma is needed
                '0,300,600,30000,30300,', # qStarts, the last comma is needed
                '0,300,600,30000,30300,', # tStarts, the last comma is needed
                ]
        writer = csv.writer(self.input_file, dialect='excel-tab')
        writer.writerow(fake_data)
        self.input_file.seek(0)

        selected_exons = parsePSL(self.input_file, 20000).next()
        self.assertEqual(len(selected_exons), 2)
        self.assertEqual(len(selected_exons[0]), 3)
        self.assertEqual(len(selected_exons[1]), 2)

    def test_remove_two_large_intron_two_exon_sets(self):
        fake_data = [
                0, # matches
                0, # misMatches
                0, # repMatches
                3, # nCount
                9, # qNumInsert
                9, # qBaseInsert
                19,# tNumInsert
                9, # tBaseInsert
                '+', # strand
                'transcript', # query sequence name
                300, # qSize
                0, # qStart
                300, # qEnd
                'chr1', # tName
                3000, # tSize
                1000, # tStart
                60100, # tEnd
                6, # blockCount
                '100,100,100,100,100,100,', # blockSize, the last comma is needed
                '0,300,600,30000,30300,60000,', # qStarts, the last comma is needed
                '0,300,600,30000,30300,60000,', # tStarts, the last comma is needed
                ]
        writer = csv.writer(self.input_file, dialect='excel-tab')
        writer.writerow(fake_data)
        self.input_file.seek(0)

        selected_exons = parsePSL(self.input_file, 20000).next()
        self.assertEqual(len(selected_exons), 2)
        self.assertEqual(len(selected_exons[0]), 3)
        self.assertEqual(len(selected_exons[1]), 2)

    def test_remove_two_large_intron_two_exon_sets(self):
        fake_data = [
                0, # matches
                0, # misMatches
                0, # repMatches
                3, # nCount
                9, # qNumInsert
                9, # qBaseInsert
                19,# tNumInsert
                9, # tBaseInsert
                '+', # strand
                'transcript', # query sequence name
                300, # qSize
                0, # qStart
                300, # qEnd
                'chr1', # tName
                3000, # tSize
                1000, # tStart
                60400, # tEnd
                7, # blockCount
                '100,100,100,100,100,100,100,', # blockSize, the last comma is needed
                '0,300,600,30000,30300,60000,60300,', # qStarts, the last comma is needed
                '0,300,600,30000,30300,60000,60300,', # tStarts, the last comma is needed
                ]
        writer = csv.writer(self.input_file, dialect='excel-tab')
        writer.writerow(fake_data)
        self.input_file.seek(0)

        selected_exons = parsePSL(self.input_file, 20000).next()
        self.assertEqual(len(selected_exons), 3)
        self.assertEqual(len(selected_exons[0]), 3)
        self.assertEqual(len(selected_exons[1]), 2)
        self.assertEqual(len(selected_exons[2]), 2)

    def test_remove_two_large_intron_flanking_one_exon(self):
        fake_data = [
                0, # matches
                0, # misMatches
                0, # repMatches
                3, # nCount
                9, # qNumInsert
                9, # qBaseInsert
                19,# tNumInsert
                9, # tBaseInsert
                '+', # strand
                'transcript', # query sequence name
                300, # qSize
                0, # qStart
                300, # qEnd
                'chr1', # tName
                3000, # tSize
                1000, # tStart
                60400, # tEnd
                6, # blockCount
                '100,100,100,100,100,100,', # blockSize, the last comma is needed
                '0,300,600,30000,60000,60300,', # qStarts, the last comma is needed
                '0,300,600,30000,60000,60300,', # tStarts, the last comma is needed
                ]
        writer = csv.writer(self.input_file, dialect='excel-tab')
        writer.writerow(fake_data)
        self.input_file.seek(0)

        selected_exons = parsePSL(self.input_file, 20000).next()
        self.assertEqual(len(selected_exons), 3)
        self.assertEqual(len(selected_exons[0]), 3)
        self.assertEqual(len(selected_exons[1]), 1)
        self.assertEqual(len(selected_exons[2]), 2)

    def tearDown(self):
        self.input_file.close()
