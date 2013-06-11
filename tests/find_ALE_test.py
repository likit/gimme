import unittest

import networkx as nx
from utils.find_ALE import Exon, find_ALE


def get_num_exons(paths):
    num_exons = []
    for path in paths:
        for exons in path:
            num_exons.append(len(exons))

    return num_exons


class TestFindALEPositive(unittest.TestCase):
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
        '''
            []-------->[]
        '''
        self.graph.add_path(['start', str(self.ex1),
                                    str(self.ex2), 'end'])
        self.transcripts = [[str(self.ex1), str(self.ex2)]]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 0)

    def test_positive_two_path_one_exon(self):
        '''
            []---------->[*]
            []------>[]
        '''

        path1 = ['start', str(self.ex1), str(self.ex2), 'end']
        path2 = ['start', str(self.ex1), str(self.ex5), 'end']
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 2)
        self.assertItemsEqual(get_num_exons(paths), [2, 2])

    def test_positive_no_ALE_with_AFE(self):
        '''
            []--------->[]--------->[]
              []------->[]--------->[]
        '''

        path1 = ['start', str(self.ex1), str(self.ex2),
                   str(self.ex3), 'end']
        path2 = ['start', str(self.ex4), str(self.ex2),
                   str(self.ex3), 'end']

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 0)

    def test_positive_single_ALE_with_AFE(self):
        '''
            []--------->[]------>[*]
              []------->[]--------->[]
        '''

        path1 = ['start', str(self.ex1), str(self.ex5),
                   str(self.ex2), 'end']
        path2 = ['start', str(self.ex1), str(self.ex5),
                   str(self.ex3), 'end']

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 2)
        self.assertItemsEqual(get_num_exons(paths), [2, 2])

    def test_positive_two_path_two_exon(self):
        '''
            []-------->[]------->[]
            []-------->[]----->[*]---->[*]
        '''
        self.ex6 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex6)] = self.ex6
        path1 = ['start', str(self.ex4),
                    str(self.ex1), str(self.ex2), 'end']
        path2 = ['start', str(self.ex4), str(self.ex1),
                    str(self.ex5), str(self.ex6), 'end']
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]

        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 2)
        self.assertItemsEqual(get_num_exons(paths), [2, 3])

    # def test_positive_two_path_multiple_exon(self):
    #     '''
    #         []---->[*]------>[*]---->[*]
    #         []-------->[]------->[]
    #     '''

    #     self.ex6 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
    #     self.ex7 = Exon('chrX', 9000, 10000, 'ex1.1', '+')
    #     self.exonsDB[str(self.ex6)] = self.ex6
    #     self.exonsDB[str(self.ex7)] = self.ex7

    #     path1 = ['start', str(self.ex1), str(self.ex3),
    #                 str(self.ex6), str(self.ex7), 'end']
    #     path2 = ['start', str(self.ex1),
    #                 str(self.ex5), str(self.ex2), 'end']
    #     self.graph.add_path(path1)
    #     self.graph.add_path(path2)
    #     self.transcripts = [path1, path2]

    #     paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

    #     self.assertEqual(len(paths[0]), 2)
    #     self.assertItemsEqual(get_num_exons(paths), [3, 4])

    def test_positive_three_path_one_exon(self):
        '''
            []---->[*]
            []-------->[**]
            []---->[]
        '''

        path1 = ['start', str(self.ex1), str(self.ex2), 'end']
        path2 = ['start', str(self.ex1), str(self.ex3), 'end']
        path3 = ['start', str(self.ex1), str(self.ex4), 'end']
        self.transcripts = [path1, path2, path3]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.graph.add_path(path3)

        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 3)
        self.assertItemsEqual(get_num_exons(paths), [2, 2, 2])

    # def test_positive_three_path_two_exon(self):
    #     '''
    #         []---->[*]---->[*]
    #         []-------->[**]---->[**]
    #         []---->[]------>[]
    #     '''
    #     self.ex6 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
    #     self.ex7 = Exon('chrX', 3500, 4500, 'ex1.1', '+')
    #     self.ex8 = Exon('chrX', 9000, 10000, 'ex1.1', '+')
    #     self.exonsDB[str(self.ex6)] = self.ex6
    #     self.exonsDB[str(self.ex7)] = self.ex7
    #     self.exonsDB[str(self.ex8)] = self.ex8

    #     path1 = ['start', str(self.ex1),
    #                 str(self.ex2), str(self.ex3), 'end']
    #     path2 = ['start', str(self.ex1),
    #                 str(self.ex5), str(self.ex6), 'end']
    #     path3 = ['start', str(self.ex1),
    #                 str(self.ex7), str(self.ex8), 'end']
    #     self.transcripts = [path1, path2, path3]

    #     self.graph.add_path(path1)
    #     self.graph.add_path(path2)
    #     self.graph.add_path(path3)

    #     paths = find_ALE(self.graph, self.exonsDB, self.transcripts)
    #     get_num_exons = [len(path) for path in paths]

    #     self.assertEqual(len(paths), 3)
    #     self.assertItemsEqual(get_num_exons, [3, 3, 3])

    def test_positive_two_path_merged_exon(self):
        '''
            []---->[*********]
            []---->[]------>[]
        '''
        self.ex8 = Exon('chrX', 2000, 5000, 'ex1.1', '+')
        self.exonsDB[str(self.ex8)] = self.ex8

        path1 = ['start', str(self.ex1),
                    str(self.ex2), str(self.ex3), 'end']
        path2 = ['start', str(self.ex1), str(self.ex8), 'end']

        self.transcripts = [path1, path2]

        self.graph.add_path(path1)
        self.graph.add_path(path2)

        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 2)
        self.assertItemsEqual(get_num_exons(paths), [3, 2])

    def test_positive_three_path_merged_exon(self):
        '''
            []-----[]---->[*********]
            []---->[]--->[]
            []-----[======]------>[]
        '''
        self.ex7 = Exon('chrX', 1000, 4000, 'ex1.1', '+')
        self.exonsDB[str(self.ex7)] = self.ex7
        self.ex8 = Exon('chrX', 2000, 5000, 'ex1.1', '+')
        self.exonsDB[str(self.ex8)] = self.ex8

        path1 = ['start', str(self.ex4), str(self.ex1),
                    str(self.ex2), 'end']
        path2 = ['start', str(self.ex4), str(self.ex1),
                    str(self.ex8), 'end']
        path3 = ['start', str(self.ex4), str(self.ex7),
                    str(self.ex3), 'end']

        self.transcripts = [path1, path2, path3]

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.graph.add_path(path3)

        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 2)
        self.assertItemsEqual(get_num_exons(paths), [2, 2])


class TestFindALENegative(unittest.TestCase):
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

    def test_negative_no_ALE_with_AFE(self):
        '''
            []<------[]<-----[]<-------[]
            []<------[]<--[*]
        '''

        path1 = (['end', str(self.ex4), str(self.ex1),
                    str(self.ex2), str(self.ex3), 'start'])
        path2 = (['end', str(self.ex4),
                    str(self.ex1), str(self.ex5), 'start'])

        path1.reverse()
        path2.reverse()

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 0)

    def test_negative_single_ALE_with_AFE(self):
        '''
               []<-------[]<-------[]
                  [*]<---[]<----------[]
        '''

        self.ex6 = Exon('chrX', 7000, 8000, 'ex1.1', '+')
        self.exonsDB[str(self.ex6)] = self.ex6

        path1 = (['end', str(self.ex4), str(self.ex2),
                    str(self.ex3), 'start'])
        path2 = (['end', str(self.ex1),
                    str(self.ex2), str(self.ex6), 'start'])

        path1.reverse()
        path2.reverse()

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 2)
        self.assertItemsEqual(get_num_exons(paths), [2, 2])

    def test_negative_one_path_one_exon(self):
        '''
            []<------[]
        '''
        path = ['end', str(self.ex1), str(self.ex2), 'start']
        path.reverse()
        self.graph.add_path(path)
        self.transcripts = [path]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 0)

    def test_negative_two_path_one_exon(self):
        '''
            []<-----------[]
                [*]<------[]
        '''
        path1 = ['end', str(self.ex1), str(self.ex2), 'start']
        path2 = ['end', str(self.ex4), str(self.ex2), 'start']

        path1.reverse()
        path2.reverse()

        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.transcripts = [path1, path2]
        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 2)
        self.assertItemsEqual(get_num_exons(paths), [2, 2])

    # def test_negative_two_path_multiple_exon(self):
    #     '''
    #         []<-----------[]-------[]
    #             [*]<--------[*]----[]
    #     '''

    #     path1 = ['end', str(self.ex1), str(self.ex2),
    #                 str(self.ex3), 'start']
    #     path2 = ['end', str(self.ex4), str(self.ex5),
    #                 str(self.ex3), 'start']
    #     path1.reverse()
    #     path2.reverse()

    #     self.graph.add_path(path1)
    #     self.graph.add_path(path2)
    #     self.transcripts = [path1, path2]

    #     paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

    #     self.assertEqual(len(paths[0]), 2)
    #     self.assertItemsEqual(get_num_exons(paths), [3, 3])

    def test_negative_three_path_one_exon(self):
        '''
            []<------------------[]
                [*]<-------------[]
                 [**]<-----------[]
        '''
        path1 = ['end', str(self.ex1), str(self.ex3), 'start']
        path2 = ['end', str(self.ex2), str(self.ex3), 'start']
        path3 = ['end', str(self.ex4), str(self.ex3), 'start']
        path1.reverse()
        path2.reverse()
        path3.reverse()

        self.transcripts = [path1, path2, path3]
        self.graph.add_path(path1)
        self.graph.add_path(path2)
        self.graph.add_path(path3)

        paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

        self.assertEqual(len(paths[0]), 3)
        self.assertItemsEqual(get_num_exons(paths), [2, 2, 2])

    # def test_negative_three_path_two_exon(self):
    #     '''
    #         []<-----------[]-------[]
    #             [*]<------[]-------[]
    #              [**]<----[]-------[]
    #     '''
    #     self.ex6 = Exon('chrX', 5500, 6500, 'ex1.1', '+')
    #     self.ex7 = Exon('chrX', 2500, 3500, 'ex1.1', '+')
    #     self.exonsDB[str(self.ex6)] = self.ex6
    #     self.exonsDB[str(self.ex7)] = self.ex7

    #     path1 = ['end', str(self.ex1), str(self.ex2), str(self.ex6), 'start']
    #     path2 = ['end', str(self.ex4), str(self.ex3), str(self.ex6), 'start']
    #     path3 = ['end', str(self.ex5), str(self.ex7), str(self.ex6), 'start']
    #     path1.reverse()
    #     path2.reverse()
    #     path3.reverse()
    #     self.transcripts = [path1, path2, path3]

    #     self.graph.add_path(path1)
    #     self.graph.add_path(path2)
    #     self.graph.add_path(path3)

    #     paths = list(find_ALE(self.graph, self.exonsDB, self.transcripts))

    #     self.assertEqual(len(paths[0]), 3)
    #     self.assertItemsEqual(get_num_exons, [3, 3, 3])
