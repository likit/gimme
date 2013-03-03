'''The script identifies mutually exclusive exons from a BED file.
Output is written in GFF format suitable for differential
exon usage analysis using MISO.

'''

import sys
import csv

import networkx as nx
from bx.intervals import IntervalTree


class Exon(object):
    def __init__(self, chrom, start, end, transcript_id, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.transcript_id = transcript_id
        self.geneID = transcript_id.split('.')[0]
        self.strand = strand

    def __str__(self):
        return "%s:%d-%d" % (self.chrom, self.start, self.end)


def parse_BED(filename):
    reader = csv.reader(open(filename), dialect='excel-tab')
    for row in reader:
        chrom = row[0]
        chrom_start = int(row[1]) + 1
        transcript_id = row[3].replace(':', '-')
        strand = row[5]
        exon_sizes = [int(s) for s in row[10].split(',')]
        exon_starts = [int(s) for s in row[11].split(',')]
        if strand == '.':
            continue
        yield (chrom,
                chrom_start,
                transcript_id,
                strand,
                exon_sizes,
                exon_starts)


def get_exon_node(infile):
    for features in parse_BED(infile):
        (chrom, chrom_start, transcript_id,
                strand, exon_sizes, exon_starts) = features
        exons = []
        for i in range(len(exon_sizes)):
            start = chrom_start + exon_starts[i]
            end = start + exon_sizes[i] - 1
            exons.append(Exon(chrom, start, end, transcript_id, strand))
        yield exons, transcript_id


def find_MXE(graph, exonsDB):
    all_exons = [exonsDB[exon] for exon in graph.nodes()]
    all_exons = sorted(all_exons, key=lambda x: x.start)
    events = []
    i = 0
    while i < len(all_exons):
        # print >> sys.stderr, str(all_exons[i]), \
        #                     graph.out_degree(str(all_exons[i]))
        if graph.out_degree(str(all_exons[i])) > 1:
            idx = i
            in_degree_found = False
            while idx < len(all_exons):
                # print >> sys.stderr, '\t', idx, str(all_exons[idx]), \
                #     graph.in_degree(str(all_exons[idx])), in_degree_found,
                if graph.in_degree(str(all_exons[idx])) > 1:
                    if not in_degree_found:
                        in_degree_found = True
                        last_idx = idx
                        if idx == len(all_exons) - 1:  # last exon
                            events.append((all_exons[i],
                                            all_exons[last_idx]))
                            # print >> sys.stderr, str(all_exons[last_idx])
                else:
                    if in_degree_found:
                        events.append((all_exons[i], all_exons[last_idx]))
                        # print >> sys.stderr, str(all_exons[last_idx])
                        in_degree_found = False

                # print >> sys.stderr, ''
                idx += 1
        i += 1
    
    # print >> sys.stderr, len(events)

    mxe_events = []
    for event in events:
        visited_path = set()
        start_exon, end_exon = event
        paths = list(nx.all_simple_paths(graph, str(start_exon), str(end_exon)))
        for path in paths:
            mxe_paths = [path]
            if '-'.join(path) not in visited_path:
                for p in paths:
                    if (p[0] == path[0] and p[-1] == path[-1]):
                        if ((p[1:-1] != path[1:-1]) and
                                p[1:-1] and path[1:-1]):
                            mxe = True  # mutually exclusive
                            for j in path[1:-1]:
                                for k in p[1:-1]:
                                    if ((j, k) in graph.edges() or
                                            (k, j) in graph.edges()):
                                        mxe = False
                            if mxe:
                                visited_path.add('-'.join(p))
                                mxe_paths.append(p)

            if len(mxe_paths) > 1:
                mxe_events.append(mxe_paths)

    # for evnt in mxe_events:
    #     for path in evnt:
    #         print >> sys.stderr, '-'.join(path)

    return mxe_events


def remove_overlaps(events, exonsDB):
    tree = IntervalTree()
    all_nodes = set()
    for path in events:
        for node in path:
            all_nodes.add(node)
            exon = exonsDB[node]
            tree.add_interval(exon)

    overlapped_exons = set()
    for node in all_nodes:
        exon = exonsDB[node]
        for overlap in tree.find(exon.start, exon.end):
            if (overlap.start != exon.start or
                    overlap.end != exon.end):
                overlapped_exons.add(node)

    new_events = []
    for path in events:
        if len(set(path).intersection(overlapped_exons)) == 0:
            new_events.append(path)
    return new_events

def write_GFF(events, exonsDB, no_events):
    all_exons = set()
    for event in events:
        for exon in event:
            all_exons.add(exonsDB[exon])
    all_exons = sorted(list(all_exons), key=lambda x: x.end)

    first_exon = all_exons[0]
    last_exon = all_exons[-1]
    mrnaid = 1
    event_no = str(no_events[first_exon.geneID])
    geneID = first_exon.geneID + '.ev' + event_no
    print "%s\tMXE\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s" % (
            first_exon.chrom, first_exon.start, last_exon.end,
            first_exon.strand, geneID, first_exon.geneID)
    for event in events:
        event_exons = sorted([exonsDB[exon] for exon in event],
                                                key=lambda x: x.end)
        first_exon = event_exons[0]
        last_exon = event_exons[-1]
        print "%s\tMXE\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.%d;Parent=%s" % (
                        first_exon.chrom, first_exon.start, last_exon.end,
                        first_exon.strand, geneID, mrnaid, geneID)
        exonid = 1
        for exon in event_exons:
            print "%s\tMXE\texon\t%d\t%d\t.\t%s\t.\tID=%s.%d.%d;Parent=%s.%d" \
                            % (exon.chrom, exon.start, exon.end,
                                exon.strand, geneID, mrnaid, exonid,
                                geneID, mrnaid)
            exonid += 1
        mrnaid += 1


def main():
    no_events = {}  # number of events in a gene
    exonsDB = {}
    infile = sys.argv[1]
    graph = nx.DiGraph()
    current_id = None
    for exons, transcript_id in get_exon_node(infile):
        new_id = transcript_id.split('.')[0]
        # print >> sys.stderr, current_id, new_id
        if not current_id:  # first gene
            for e in exons:
                exonsDB[str(e)] = e
            graph.add_path([str(e) for e in exons])
            current_id = new_id
            no_events[current_id] = 0
        else:
            if new_id != current_id:
                for events in find_MXE(graph, exonsDB):
                    events = remove_overlaps(events, exonsDB)
                    if len(events) > 1:
                        no_events[current_id] += 1
                        write_GFF(events, exonsDB, no_events)

                graph = nx.DiGraph()
                exonsDB = {}
                current_id = new_id
                no_events[current_id] = 0

            for e in exons:
                exonsDB[str(e)] = e
            graph.add_path([str(e) for e in exons])

    for events in find_MXE(graph, exonsDB):
        events = remove_overlaps(events, exonsDB)
        if len(events) > 1:
            no_events[current_id] += 1
            write_GFF(events, exonsDB, no_events)

if __name__ == '__main__':
    main()
