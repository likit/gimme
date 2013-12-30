'''
    Gimme: a transcripts assembler based on alignments.

    Copyright (C) 2012 Michigan State University

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Likit Preeyanon, preeyano@msu.edu
'''
#!/usr/bin/env python

import sys
import csv
import argparse

from sys import stderr, stdout

import networkx as nx

#from matplotlib import pyplot as plt
from utils import pslparser, get_min_isoforms, split_strand
from bx.intervals.intersection import Interval, IntervalTree
from pygr import seqdb


gap_size = 50  # a minimum intron size (bp)
max_intron = 300000  # a maximum intron size (bp)
min_utr = 100  # a minimum UTR size (bp)
min_transcript_len = 300  # a minimum length for multiple exon transcript(bp)
min_single_exon_len = 500  # a minimum length for a single exon(bp)
max_isoforms = 20   # minimal isoforms will be searched
                    #if the number of isoforms exceed this number
VERSION = '0.98'


class ExonObj:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.terminal = None
        self.next_exons = set()
        self.introns = set()
        self.single = False
        self.remove = False

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)

    def get_size(self):
        return self.end - self.start + 1


class AlignmentDB(object):
    def __init__(self):
        self.exon_db = {}  # store all exon objects
        self.intron_db = {}  # store all intron objects
        self.single_exons_db = {}  # store all single exon objects
        self.single_exons_intervals = {}  # store intersecter objects for
                                          # single exons


def parse_bed(bed_file):
    '''Reads alignments from BED format and creates
    exon objects from a transcript.

    '''
    reader = csv.reader(bed_file, dialect='excel-tab')
    for row in reader:
        exons = []
        chrom = row[0]
        chrom_start = int(row[1])

        exon_sizes = [int(s) for s in row[10].split(',')]
        exon_starts = [chrom_start + int(s) for s in row[11].split(',')]

        for i in range(len(exon_starts)):
            exon_start = exon_starts[i]
            exon_end = exon_start + exon_sizes[i]

            exon = ExonObj(chrom, exon_start, exon_end)
            exons.append(exon)

        exons = delete_gap(exons, gap_size)
        yield exons


def parse_psl(psl_file):
    '''Reads alignments from PSL format and creates
    exon objects from each transcript.

    '''
    for pslobj in pslparser.read(psl_file):
        exons = []

        for i in range(len(pslobj.tStarts)):
            exon_start = pslobj.tStarts[i]
            exon_end = exon_start + pslobj.blockSizes[i]

            exon = ExonObj(pslobj.tName, exon_start, exon_end)
            exons.append(exon)

        exons = delete_gap(exons, gap_size)
        yield exons


def remove_large_intron(exons, max_intron=1e6):
    '''Returns groups of exons split by introns longer than
    max_intron or return a list containing the original exons.

    Exons is a list of exons sorted by start position.

    '''
    if max_intron < 0:
        return [exons]  # disable

    all_exon_groups = []
    subgroup = []
    for i in range(len(exons)):
        curr_exon = exons[i]
        try:
            next_exon = exons[i + 1]
        except IndexError:  # end of list
            subgroup.append(curr_exon)
            all_exon_groups.append(subgroup[:])
            break
        else:
            intron_start = curr_exon.end + 1
            intron_end = next_exon.start - 1

            if intron_end - intron_start > max_intron:
                subgroup.append(curr_exon)
                all_exon_groups.append(subgroup[:])

                subgroup = []  # start new subgroup
            else:
                subgroup.append(curr_exon)

    return all_exon_groups


def add_intron(exons, align_db, clusters, cluster_no):
    '''Get introns from a set of exons.

    Intron object is created for each intron with all connected exons
    as an attribute.
    '''

    introns = []
    existing_clusters = set()

    for i in range(len(exons)):
        curr_exon = exons[i]
        try:
            next_exon = exons[i + 1]
        except IndexError:
            pass
        else:
            intron_start = curr_exon.end + 1
            intron_end = next_exon.start - 1

            curr_exon.next_exons.add(str(next_exon))

            intron_name = '%s:%d-%d' % (curr_exon.chrom,
                                            intron_start,
                                            intron_end)
            intron = nx.DiGraph(name=intron_name, cluster=None)

            try:
                intron_ = align_db.intron_db[intron.graph['name']]
            except KeyError:
                align_db.intron_db[intron.graph['name']] = intron
                intron.add_edge(str(curr_exon), str(next_exon))
                introns.append(intron)

                curr_exon.introns.add(intron.graph['name'])
                next_exon.introns.add(intron.graph['name'])
            else:
                intron_.add_edge(str(curr_exon), str(next_exon))
                introns.append(intron_)
                existing_clusters.add(intron_.graph['cluster'])

                curr_exon.introns.add(intron_.graph['name'])
                next_exon.introns.add(intron_.graph['name'])

    if introns:
        cluster_no += 1  # create new cluster index
        if not existing_clusters:
            cluster = nx.DiGraph()
            if len(introns) > 1:
                cluster.add_path([i.graph['name'] for i in introns])
                for intron in introns:
                    intron.graph['cluster'] = cluster_no
            else:
                cluster.add_node(introns[0].graph['name'])
                introns[0].graph['cluster'] = cluster_no
        else:
            cluster = nx.DiGraph(exons=set())
            for cl in existing_clusters:
                cluster.add_edges_from(clusters[cl].edges())
                clusters.pop(cl)

            for intron in cluster.nodes():
                align_db.intron_db[intron].graph['cluster'] = cluster_no

            if len(introns) > 1:
                cluster.add_path([i.graph['name'] for i in introns])

                for intron in introns:
                    intron.graph['cluster'] = cluster_no
            else:
                cluster.add_node(introns[0].graph['name'])
                introns[0].graph['cluster'] = cluster_no

        clusters[cluster_no] = cluster

    return cluster_no


def collapse_exon(g, align_db):
    '''Merge overlapped exons together.

    An exon gets extended when they are merged with a larger exon.

    A smaller exon is then removed from the graph.

    '''

    exons = [align_db.exon_db[e] for e in g.nodes()]
    sorted_exons = sorted(exons, key=lambda x: (x.end, x.start))
    chromosome = exons[0].chrom

    i = 0
    curr_exon = sorted_exons[i]
    while i <= len(sorted_exons):
        try:
            next_exon = sorted_exons[i + 1]
        except IndexError:
            pass
        else:
            if curr_exon.end == next_exon.end:
                if next_exon.terminal == 1:  # left terminal
                    g.add_edges_from([(str(curr_exon), n)\
                            for n in g.successors(str(next_exon))])
                    g.remove_node(str(next_exon))
                    if curr_exon.terminal == 2:
                        curr_exon.terminal = None
                else:
                    if (curr_exon.terminal == 1 and
                            next_exon.start - curr_exon.start <= min_utr):
                        g.add_edges_from([(str(next_exon), n) for n in
                                            g.successors(str(curr_exon))])
                        g.remove_node(str(curr_exon))
                    curr_exon = next_exon
            else:
                curr_exon = next_exon
        i += 1

    i = 0
    exons = [align_db.exon_db[e] for e in g.nodes()]
    sorted_exons = sorted(exons, key=lambda x: (x.start, x.end))
    curr_exon = sorted_exons[0]
    while i <= len(sorted_exons):
        try:
            next_exon = sorted_exons[i + 1]
        except IndexError:
            pass
        else:
            if curr_exon.start == next_exon.start:
                if curr_exon.terminal == 2:
                    g.add_edges_from([(n, str(next_exon))\
                            for n in g.predecessors(str(curr_exon))])
                    g.remove_node(str(curr_exon))
                    curr_exon = next_exon
                else:
                    if next_exon.terminal == 2:
                        if next_exon.end - curr_exon.end <= min_utr:
                            g.add_edges_from([(n, str(curr_exon))\
                                    for n in g.predecessors(str(next_exon))])
                            g.remove_node(str(next_exon))
                        else:
                            curr_exon = next_exon
                    else:
                        curr_exon = next_exon
            else:
                curr_exon = next_exon
        i += 1
    try:
        '''If there are single exons in this chromosome,
        remove or extend them according to how they overlap with
        exons in the cluster.

        '''
        singles = align_db.single_exons_intervals[chromosome]
    except KeyError:
        pass
    else:
        for node in g.nodes():
            exon = align_db.exon_db[node]
            remove_redundant_exon(exon, singles, set())


def remove_redundant_exon(exon, singles, unmergables):
    '''Recursively remove exons that are subset of a given exon.'''

    overlaps = [o for o in singles.find(exon.start, exon.end) \
                                if not o.value['exon'].remove and \
                                str(o.value['exon']) not in unmergables]
    if not overlaps:
        return
    for o in overlaps:
        if o.start >= exon.start and o.end <= exon.end:
            o.value['exon'].remove = True  # mark the exon as removed
        elif o.start >= exon.start and o.end > exon.end:
            if o.end - exon.end < min_utr:
                o.value['exon'].remove = True
            else:
                unmergables.add(str(o.value['exon']))
        elif o.start < exon.start and o.end <= exon.end:
            if exon.start - o.start < min_utr:
                o.value['exon'].remove = True
            else:
                unmergables.add(str(o.value['exon']))
        elif o.start < exon.start and o.end > exon.end:
            if (exon.start - o.start) + (o.end - exon.end) < min_utr:
                o.value['exon'].remove = True
            else:
                unmergables.add(str(o.value['exon']))
        else:
            unmergables.add(str(o.value['exon']))

    remove_redundant_exon(exon, singles, unmergables)


def delete_gap(exons, gap_size=0):
    '''Alignments may contain small gaps from indels and etc.

    The program fills up gaps to obtain a complete exon.

    A maximum size of a gap can be adjusted by assigning a new
    value to gap_size parameter on a command line.

    '''

    i = 0
    new_exons = []
    curr_exon = exons[i]

    while True:
        try:
            next_exon = exons[i + 1]
        except IndexError:
            break
        else:
            if next_exon.start - curr_exon.end <= gap_size:
                curr_exon.end = next_exon.end
            else:
                new_exons.append(curr_exon)
                curr_exon = next_exon
        i += 1

    new_exons.append(curr_exon)
    return new_exons


def add_exon(align_db, exons):
    '''1.Change a terminal attribute of a leftmost exon
    and a rightmost exon to 1 and 2 respectively.

    A terminal attribute has a value 'None' by default.

    2.Add exons to the exon database (db).

    '''
    exons[0].terminal = 1  # left end
    exons[-1].terminal = 2  # right end

    for exon in exons:
        try:
            exon_ = align_db.exon_db[str(exon)]
        except KeyError:
            align_db.exon_db[str(exon)] = exon
        else:
            if ((exon.terminal and exon_.terminal) and
                        exon.terminal != exon_.terminal):
                exon.terminal = None

            if not exon.terminal and exon_.terminal:
                exon_.terminal = None


def merge_cluster(align_db):
    '''Connect introns from the same gene together.'''

    big_cluster = nx.Graph()
    paths = []
    for exon in align_db.exon_db.itervalues():
        pth = []
        for intron in exon.introns:
            cluster = align_db.intron_db[intron].graph['cluster']
            pth.append(cluster)
            # exon.clusters.add(cluster)
        paths.append(pth)

    for pth in paths:
        if len(pth) > 2:
            big_cluster.add_path(pth)
        elif len(pth) == 2:
            big_cluster.add_edges_from([pth])
        elif len(pth) == 1:
            big_cluster.add_nodes_from(pth)
        else:
            pass

    return big_cluster


def print_bed(align_db, transcript, strand, gene_id, tran_id):
    '''Print a splice graph in BED format.'''

    exons = [align_db.exon_db[e] for e in transcript]

    chrom_start = exons[0].start
    chrom_end = exons[-1].end
    chrom = exons[0].chrom

    block_starts = ','.join([str(exon.start - chrom_start) for exon in exons])
    block_sizes = ','.join([str(exon.end - exon.start) for exon in exons])

    name = '%s:%d.%d' % (chrom, gene_id, tran_id)
    score = 1000
    item_RGB = '0,0,0'
    thick_start = chrom_start
    thick_end = chrom_end
    block_count = len(exons)

    writer = csv.writer(stdout, dialect='excel-tab')
    writer.writerow((chrom,
                    chrom_start,
                    chrom_end,
                    name,
                    score,
                    strand,
                    thick_start,
                    thick_end,
                    item_RGB,
                    block_count,
                    block_sizes,
                    block_starts))


def print_bed_single(exon, gene_id, tran_id):
    '''Print a splice graph in BED format.'''

    chrom_start = exon.start
    chrom_end = exon.end
    chrom = exon.chrom

    block_starts = ','.join([str(exon.start - chrom_start)])
    block_sizes = ','.join([str(exon.end - exon.start)])

    name = '%s:%d.%d' % (chrom, gene_id, tran_id)
    score = 1000
    item_RGB = '0,0,0'
    thick_start = chrom_start
    thick_end = chrom_end
    strand = '+'
    block_count = 1

    writer = csv.writer(stdout, dialect='excel-tab')
    writer.writerow((chrom,
                    chrom_start,
                    chrom_end,
                    name,
                    score,
                    strand,
                    thick_start,
                    thick_end,
                    item_RGB,
                    block_count,
                    block_sizes,
                    block_starts))


def build_gene_model(genome,
                        align_db,
                        clusters,
                        big_cluster,
                        find_max,
                        min_transcript_len=0,
                        max_isoforms=1e6,
                    ):

    '''Build and print out gene models.'''

    visited_clusters = set()
    transcripts_num = 0
    gene_id = 0
    excluded = 0
    two_exon_trns = set()

    def check_criteria(transcript, two_exon_trns):
        '''Return True or False whether a transcript pass or
        fail the criteria.

        '''
        transcript_length = sum([align_db.exon_db[e].get_size() \
                                                for e in transcript])

        if transcript_length <= min_transcript_len:
            return False  # fail
        else:
            if len(transcript) == 2:
                trns = ','.join(transcript)
                if trns in two_exon_trns:
                    return False  # fail
                else:
                    two_exon_trns.add(trns)
                    return True  # pass
            else:
                return True

    def exon_to_exonobj(exon):
        '''Returns an exon objects from a given exon coordinate.'''
        chrom, coord = exon.split(':')
        start, end = coord.split('-')
        return ExonObj(chrom, int(start), int(end))

    for cl_num, cl in enumerate(big_cluster.nodes(), start=1):
        if cl not in visited_clusters:
            g = nx.DiGraph()
            for intron in clusters[cl].nodes():
                g.add_edges_from(align_db.intron_db[intron].edges())

            visited_clusters.add(cl)

            for neighbor in nx.dfs_tree(big_cluster, cl):
                neighbor_cluster = clusters[neighbor]
                for intron in neighbor_cluster.nodes():
                    g.add_edges_from(align_db.intron_db[intron].edges())

                visited_clusters.add(neighbor)
            # # nx.draw_spring(nx.algorithms.dfs_tree(g))
            # nx.draw_spring(g)
            # plt.show()
            # for node in g.nodes():
            #     print node, g[node]
            # raise SystemExit
            collapse_exon(g, align_db)
            for g in split_strand.split(g, genome):
                if g.nodes():
                    subalign_db = AlignmentDB()
                    for edge in g.edges():
                        exon1 = exon_to_exonobj(edge[0])
                        exon2 = exon_to_exonobj(edge[1])
                        add_exon(subalign_db, [exon1, exon2])
                    collapse_exon(g, subalign_db)

                    trans_id = 0
                    gene_id += 1
                    strand = g.graph['strand']
                    for node in g.nodes():
                        if not g.predecessors(node):
                            g.add_edge('Start', node)
                        if not g.successors(node):
                            g.add_edge(node, 'End')

                    max_paths = [path for path in \
                                    nx.all_simple_paths(g, 'Start', 'End')]

                    if find_max:
                        '''Report all maximum isoforms.'''

                        for transcript in max_paths:
                            transcript = transcript[1:-1]
                            if check_criteria(transcript, two_exon_trns):
                                transcripts_num += 1
                                trans_id += 1
                                print_bed(align_db,
                                            transcript,
                                            strand,
                                            gene_id,
                                            trans_id)
                            else:
                                excluded += 1
                    else:
                        '''Report minimal isoforms if maximum isoforms exceeds
                        max_isoforms.

                        '''
                        if len(max_paths) > max_isoforms:
                            for transcript in \
                                    get_min_isoforms.get_min_paths(g, False):
                                if check_criteria(transcript, two_exon_trns):
                                    transcripts_num += 1
                                    trans_id += 1
                                    print_bed(align_db,
                                                transcript,
                                                strand,
                                                gene_id,
                                                trans_id)
                                else:
                                    excluded += 1
                        else:
                            for transcript in max_paths:
                                transcript = transcript[1:-1]
                                if check_criteria(transcript, two_exon_trns):
                                    transcripts_num += 1
                                    trans_id += 1
                                    print_bed(align_db,
                                                transcript,
                                                strand,
                                                gene_id,
                                                trans_id)
                                else:
                                    excluded += 1

        print >> stderr, '\r  |--Multi-exon\t\t%d genes, %d isoforms ' % \
                                            (gene_id, transcripts_num),

    return gene_id, transcripts_num, excluded


def merge_exon(align_db):
    '''Return merged exons from exons overlapped to each other.'''

    exons = {}
    for chrom in align_db.single_exons_db:
        exons[chrom] = []
        for exn in align_db.single_exons_db[chrom]:
            if not exn.remove:
                exons[chrom].append(exn)
        exons[chrom] = sorted(exons[chrom], key=lambda x: x.start)

    new_exons = {}

    for chrom in exons:
        i = 0
        new_exons[chrom] = []
        curr_exon = exons[chrom][i]
        while i < len(exons[chrom]):
            try:
                next_exon = exons[chrom][i + 1]
            except IndexError:
                new_exons[chrom].append(curr_exon)
                break
            else:
                if next_exon.start <= curr_exon.end:
                    if next_exon.end > curr_exon.end:
                        next_exon.start = curr_exon.start
                        curr_exon = next_exon
                else:
                    new_exons[chrom].append(curr_exon)
                    curr_exon = next_exon
            i += 1

    return new_exons


def detect_format(input_file):
    '''Returns a file format detected from input file.'''

    fp = open(input_file)
    cols = fp.readline().split()
    fp.close()

    if len(cols) == 21:
        if int(cols[11]) <= int(cols[12]) and cols[8] in ['+', '.', '-']:
            return 'PSL'
    elif len(cols) == 12:
        if int(cols[1]) <= int(cols[2]) and cols[5] in ['+', '.', '-']:
            return 'BED'
    else:
        return None


def main(input_files):
    print >> stderr, 'Gimme : Alignment-based assembler'
    print >> stderr, 'Version : %s' % (VERSION)
    print >> stderr, 'Source code : https://github.com/ged-lab/gimme.git\n'
    print >> stderr, 'Building a sequence DB...'
    genome = seqdb.SequenceFileDB(args.reference)

    if args.debug:
        print >> stderr, 'DEBBUG MODE\t' + \
                'Use this mode for debugging only!\n'

    print >> stderr, '[Run...]'

    cluster_no = 0
    clusters = {}
    align_db = AlignmentDB()

    '''======Detect input format======'''
    for input_file in input_files:
        input_format = detect_format(input_file)
        if input_format == 'PSL':
            parse = parse_psl
        elif input_format == 'BED':
            parse = parse_bed
        else:
            print >> stderr, 'ERROR: Unrecognized input format. ' + \
                    'Use utils/gff2bed.py to convert GFF to BED.'
            raise SystemExit

        '''====Parse alignments and build exon objects===='''
        print >> stderr, 'Input\t\t\t%s' % input_file
        for n, exons in enumerate(parse(open(input_file)), start=1):
            for group in remove_large_intron(exons, max_intron):
                if len(group) > 1:
                    add_exon(align_db, group)  # add exons to exon db
                    cluster_no = add_intron(group, align_db,
                                                clusters, cluster_no)
                else:
                    exon = group[0]  # add a lone exon to single exon db
                    if exon.chrom not in align_db.single_exons_db:
                        align_db.single_exons_db[exon.chrom] = [exon]
                    else:
                        align_db.single_exons_db[exon.chrom].append(exon)

            if n % 100 == 0:
                print >> stderr, '\r  |--Parsing\t\t%d alignments' % n,
        print >> stderr, '\r  |--Parsing\t\t%d alignments' % n

    '''====Merge overlapped single exons===='''
    merged_single_exons = merge_exon(align_db)

    '''====Build intervals from single exons.===='''
    for chrom in merged_single_exons:
        align_db.single_exons_intervals[chrom] = IntervalTree()
        for exon in merged_single_exons[chrom]:
            interval = Interval(exon.start, exon.end, value={'exon': exon})
            align_db.single_exons_intervals[chrom].insert_interval(interval)

    '''====Connect introns from the same gene to each other===='''
    big_cluster = merge_cluster(align_db)

    '''====Build gene models===='''
    print >> stderr, 'Constructing'
    return_items = build_gene_model(genome,
                                        align_db,
                                        clusters,
                                        big_cluster,
                                        args.max,
                                        min_transcript_len,
                                        max_isoforms,
                                    )

    print >> stderr, ''
    gene_id, transcripts_num, excluded = return_items

    single_exon_gene_num = 0
    for chrom in merged_single_exons:
        for exon in merged_single_exons[chrom]:
            if (exon.get_size() > min_single_exon_len
                                    and not exon.remove):
                gene_id += 1
                transcripts_num += 1
                single_exon_gene_num += 1
                print_bed_single(exon, gene_id, 1)
                print >> stderr, '\r  |--Single-exon\t%d genes' % \
                                                    single_exon_gene_num,
            else:
                excluded += 1

    '''====Print out summary report to standard error===='''
    print >> stderr, '\n[Done]'
    if gene_id > 0:
        print >> stderr, \
            '\nTotal multi-exon gene = %d gene(s) / %d isoform(s)' % \
                                                (gene_id, transcripts_num)
        print >> stderr, \
            'Total single-exon gene = %d gene(s)' % single_exon_gene_num
    else:
        print >> stderr, 'No gene models built.',
    if excluded > 0 and args.debug:
        print >> stderr, '(%d transcripts do not pass criteria.)' % excluded
    else:
        print >> stderr, ''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='gimme.py')
    parser.add_argument('--min_utr', type=int, metavar='int',
            default=min_utr,
            help='a cutoff size of alternative UTRs (bp)' +
                    ' (default: %(default)s)')
    parser.add_argument('--gap_size', type=int, metavar='int',
            default=gap_size,
            help='the maximum gap size (bp) (default: %(default)s)')
    parser.add_argument('--max_intron', type=int, metavar='int',
            default=max_intron,
            help='the maximum intron size (bp) (default: %(default)s)')
    parser.add_argument('--max_isoforms', type=int, metavar='int',
            default=max_isoforms,
            help='the maximum number of isoforms reported ' +
            'without -x option (default: %(default)s)')
    parser.add_argument('--min_transcript_len', type=int,
            metavar='int', default=min_transcript_len,
            help='the minimum size of transcript (bp)' +
                    '(default: %(default)s)')
    parser.add_argument('--min_single_exon_len', type=int,
            metavar='int', default=min_single_exon_len,
            help='the minimum size of a transcript with a single exon (bp)' +
                    '(default: %(default)s)')
    parser.add_argument('-x', '--max', action='store_true',
            help='report all putative isoforms')
    parser.add_argument('--debug', action='store_true',
            help='reset parameters (for debugging purpose only)')
    parser.add_argument('input', type=str, nargs='+',
            help='input file(s) in PSL/BED format')
    parser.add_argument('-v', '--version', action='version',
            version='%(prog)s version ' + VERSION)
    parser.add_argument('-r','--reference', type=str,
            help='a reference genome in FASTA format')

    args = parser.parse_args()
    if not args.reference:
        print >> sys.stderr, "A reference file is required."
        sys.exit()

    if args.debug:
        '''Parameters are set to retain all splice junctions for
        debugging.

        '''
        gap_size = 0
        max_intron = -1
        min_utr = 0
        min_transcript_len = 1
        min_single_exon_len = 1
        args.max = True
    else:
        if args.min_utr <= 0:
            raise ValueError('Invalid UTRs size (<=0)')
        elif args.min_utr != min_utr:
            min_utr = args.min_utr
            print >> sys.stderr, 'User defined min_utr = %d' % min_utr

        if args.gap_size < 0:
            raise ValueError('Invalid intron size (<0)')
        elif args.gap_size != gap_size:
            gap_size = args.gap_size
            print >> sys.stderr, 'User defined gap_size = %d' % gap_size

        if args.max_intron <= 0:
            raise ValueError('Invalid intron size (<=0)')
        elif args.max_intron != max_intron:
            max_intron = args.max_intron
            print >> sys.stderr, 'User defined max_intron = %d' % max_intron

        if args.max_isoforms <= 0:
            raise ValueError('Invalid number of isoforms (<=0)')
        elif args.max_isoforms != max_isoforms:
            max_isoforms = args.max_isoforms
            print >> sys.stderr, \
                    'User defined max_isoforms = %d' % max_isoforms

        if args.min_transcript_len <= 0:
            raise ValueError('Invalid transcript size (<=0)')
        elif args.min_transcript_len != min_transcript_len:
            min_transcript_len = args.min_transcript_len
            print >> sys.stderr, 'User defined min_transcript_len = %d' % \
                                                        min_transcript_len
        if args.min_single_exon_len <= 0:
            raise ValueError('Invalid transcript size (<=0)')
        elif args.min_single_exon_len != min_single_exon_len:
            min_single_exon_len = args.min_single_exon_len
            print >> sys.stderr, 'User defined min_single_exon_len = %d' % \
                                                        min_single_exon_len
    if args.input:
        main(args.input)
