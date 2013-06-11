import unittest

import networkx as nx
from utils.find_MXE import find_MXE, Exon, remove_overlaps

class TestRemoveOverlaps(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.graph = nx.DiGraph()
        self.ex1 = Exon('chrX', 100, 200, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 300, 400, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 500, 600, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 700, 800, 'ex1.1', '+')
        self.ex5 = Exon('chrX', 900, 1000, 'ex1.1', '+')
        self.ex6 = Exon('chrX', 1100, 1200, 'ex1.1', '+')
        self.ex7 = Exon('chrX', 350, 400, 'ex1.1', '+')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4
        self.exonsDB[str(self.ex5)] = self.ex5
        self.exonsDB[str(self.ex6)] = self.ex6
        self.exonsDB[str(self.ex7)] = self.ex7

    def test_no_events(self):
        self.events = []
        self.new_events = remove_overlaps(self.events, self.exonsDB)

        self.assertEqual(len(self.new_events), 0)

    def test_no_event_one_overlap(self):
        self.path1 = [str(self.ex1), str(self.ex2), str(self.ex3)]
        self.path2 = [str(self.ex1), str(self.ex7), str(self.ex3)]
        self.events = [self.path1, self.path2]
        self.new_events = remove_overlaps(self.events, self.exonsDB)

        self.assertEqual(len(self.new_events), 0)

    def test_one_event_one_overlap(self):
        self.path1 = [str(self.ex1), str(self.ex2), str(self.ex4)]
        self.path2 = [str(self.ex1), str(self.ex7), str(self.ex4)]
        self.path3 = [str(self.ex1), str(self.ex3), str(self.ex4)]
        self.events = [self.path1, self.path2, self.path3]
        self.new_events = remove_overlaps(self.events, self.exonsDB)

        '''new_events with one path will be not included in the output'''

        self.assertEqual(len(self.new_events) - 1, 0)

class TestFindMXE(unittest.TestCase):
    def setUp(self):
        self.exonsDB = {}
        self.graph = nx.DiGraph()
        self.ex1 = Exon('chrX', 100, 200, 'ex1.1', '+')
        self.ex2 = Exon('chrX', 300, 400, 'ex1.1', '+')
        self.ex3 = Exon('chrX', 500, 600, 'ex1.1', '+')
        self.ex4 = Exon('chrX', 700, 800, 'ex1.1', '+')
        self.ex5 = Exon('chrX', 900, 1000, 'ex1.1', '+')
        self.ex6 = Exon('chrX', 1100, 1200, 'ex1.1', '+')
        self.ex7 = Exon('chrX', 1300, 1400, 'ex1.1', '+')
        self.exonsDB[str(self.ex1)] = self.ex1
        self.exonsDB[str(self.ex2)] = self.ex2
        self.exonsDB[str(self.ex3)] = self.ex3
        self.exonsDB[str(self.ex4)] = self.ex4
        self.exonsDB[str(self.ex5)] = self.ex5
        self.exonsDB[str(self.ex6)] = self.ex6
        self.exonsDB[str(self.ex7)] = self.ex7

    def test_no_event(self):
        pass
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex3),
                                str(self.ex4), str(self.ex5)])

        self.events = list(find_MXE(self.graph, self.exonsDB))

        self.assertEqual(len(self.events), 0)

    def test_no_event_one_skipped_exon(self):
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex3), str(self.ex4), str(self.ex5)])
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex4), str(self.ex5)])
        self.graph.add_path([str(self.ex1), str(self.ex3),
                                str(self.ex4), str(self.ex5)])

        self.events = list(find_MXE(self.graph, self.exonsDB))

        self.assertEqual(len(self.events), 0)

    def test_no_event_with_one_two_skipped_exon(self):
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex3), str(self.ex4), str(self.ex5)])
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex5)])

        self.events = list(find_MXE(self.graph, self.exonsDB))

        self.assertEqual(len(self.events), 0)

    def test_one_event_one_mx_exon(self):
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex3), str(self.ex6),
                                str(self.ex7)])
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex4), str(self.ex6),
                                str(self.ex7)])

        self.events = find_MXE(self.graph, self.exonsDB)
        self.assertEqual(len(self.events), 1)
        self.assertEqual(len(self.events[0]), 2)

        self.assertEqual(self.exonsDB[self.events[0][0][0]].start, 300)
        self.assertEqual(self.exonsDB[self.events[0][0][0]].end, 400)
        self.assertEqual(self.exonsDB[self.events[0][0][-1]].start, 1100)
        self.assertEqual(self.exonsDB[self.events[0][0][-1]].end, 1200)

    def test_one_event_two_mx_exon(self):
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex3), str(self.ex6),
                                str(self.ex7)])
        self.graph.add_path([str(self.ex1), str(self.ex4),
                                str(self.ex5), str(self.ex6),
                                str(self.ex7)])

        self.events = find_MXE(self.graph, self.exonsDB)
        self.assertEqual(len(self.events), 1)
        self.assertEqual(len(self.events[0]), 2)

        self.assertEqual(self.exonsDB[self.events[0][0][0]].start, 100)
        self.assertEqual(self.exonsDB[self.events[0][0][0]].end, 200)
        self.assertEqual(self.exonsDB[self.events[0][0][-1]].start, 1100)
        self.assertEqual(self.exonsDB[self.events[0][0][-1]].end, 1200)

    def test_one_event_two_mx_exon_with_three_isoforms(self):
        self.graph.add_path([str(self.ex1), str(self.ex2),
                                str(self.ex3), str(self.ex6), str(self.ex7)])
        self.graph.add_path([str(self.ex1), str(self.ex4),
                                str(self.ex5), str(self.ex6), str(self.ex7)])
        self.graph.add_path([str(self.ex1), str(self.ex2), str(self.ex4),
                                str(self.ex5), str(self.ex6), str(self.ex7)])

        self.events = find_MXE(self.graph, self.exonsDB)
        self.assertEqual(len(self.events), 1)
        self.assertEqual(len(self.events[0]), 2)

        self.assertEqual(self.exonsDB[self.events[0][0][0]].start, 300)
        self.assertEqual(self.exonsDB[self.events[0][0][0]].end, 400)
        self.assertEqual(self.exonsDB[self.events[0][0][-1]].start, 1100)
        self.assertEqual(self.exonsDB[self.events[0][0][-1]].end, 1200)
