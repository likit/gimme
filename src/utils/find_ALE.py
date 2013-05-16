import csv
import sys

import networkx as nx


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


def find_ALE(graph, exonsDB, transcripts):
    degrees = set()
    for path in list(nx.all_simple_paths(graph, 'start', 'end')):
        # print >> sys.stderr, path
        for exon in path:
            if exon == 'start' or exon == 'end':
                continue
            if graph.out_degree(exon) > 1:
                degrees.add(exon)
                break

    # print >> sys.stderr, degrees

    if not degrees:
        return []

    degrees = sorted(degrees, key=lambda x:exonsDB[x].start)

    common_exon = degrees[-1]

    paths = []
    ALE = set()
    for tranx in transcripts:
        start = common_exon
        end = tranx[-2]
        # print >> sys.stderr, tranx
        # print >> sys.stderr, 'start = ', start, 'end= ', end
        try:
            for path in list(nx.all_simple_paths(graph, start, end)):
                if path[-1] not in ALE:
                    ALE.add(path[-1])
                    paths.append(path)
        except IndexError:
            pass
    # print >> sys.stderr, 'ALE =', ALE

    if len(paths) > 1:
        # print >> sys.stderr, 'PATHS=', paths
        return paths
    else:
        return []


def write_GFF(events, exonsDB, no_events):
    all_exons = set()
    # print events
    for event in events:
        for exon in event:
            all_exons.add(exonsDB[exon])
    all_exons = sorted(list(all_exons), key=lambda x: x.end)

    first_exon = all_exons[0]
    last_exon = all_exons[-1]
    mrnaid = 1
    event_no = str(no_events[first_exon.geneID])
    geneID = first_exon.geneID + '.ev' + event_no
    print "%s\tSE\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s" % (
            first_exon.chrom, first_exon.start, last_exon.end,
            first_exon.strand, geneID, first_exon.geneID)
    for event in events:
        event_exons = sorted([exonsDB[exon] for exon in event],
                                                key=lambda x: x.end)
        first_exon = event_exons[0]
        last_exon = event_exons[-1]
        print "%s\tSE\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.%d;Parent=%s" % (
                        first_exon.chrom, first_exon.start, last_exon.end,
                        first_exon.strand, geneID, mrnaid, geneID)
        exonid = 1
        for exon in event_exons:
            print "%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s.%d.%d;Parent=%s.%d" \
                            % (exon.chrom, exon.start, exon.end,
                                exon.strand, geneID, mrnaid, exonid,
                                geneID, mrnaid)
            exonid += 1
        mrnaid += 1

    # first_exon = all_exons[0]
    # last_exon = all_exons[-1]
    # print "%s\tSE\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.%d;Parent=%s" % (
    #                 first_exon.chrom, first_exon.start, last_exon.end,
    #                 first_exon.strand, geneID, mrnaid, geneID)
    # print "%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s.%d.%d;Parent=%s.%d" % \
    #                 (exon.chrom, first_exon.start, first_exon.end,
    #                     first_exon.strand, geneID, mrnaid, 1,
    #                     geneID, mrnaid)
    # print "%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s.%d.%d;Parent=%s.%d" % \
    #                 (exon.chrom, last_exon.start, last_exon.end,
    #                     last_exon.strand, geneID, mrnaid, 2,
    #                     geneID, mrnaid)
def add_exons(exonsDB, exons, graph, transcripts):
        for e in exons:
            exonsDB[str(e)] = e

        strand = exons[0].strand

        if strand == '+':
            path = ['start']
        else:
            path = ['end']

        for e in exons:
            path.append(str(e))

        if strand == '+':
            path.append('end')
        else:
            path.append('start')
            path.reverse()

        graph.add_path(path)
        transcripts.append(path)


def main():
    no_events = {}  # number of events in a gene
    exonsDB = {}
    infile = sys.argv[1]
    graph = nx.DiGraph()
    current_id = None
    transcripts = []
    for exons, transcript_id in get_exon_node(infile):
        if len(exons) == 1:
            continue
        new_id = transcript_id.split('.')[0]
        if not current_id:  # first gene
            add_exons(exonsDB, exons, graph, transcripts)

            current_id = new_id
            no_events[current_id] = 0
        else:
            if new_id != current_id:
                if len(transcripts) > 1:
                    events = find_ALE(graph, exonsDB, transcripts)
                    if events:
                        no_events[current_id] += 1
                        write_GFF(events, exonsDB, no_events)

                graph = nx.DiGraph()
                exonsDB = {}
                current_id = new_id
                no_events[current_id] = 0
                transcripts = []

            add_exons(exonsDB, exons, graph, transcripts)

    if len(transcripts) > 1:
        events = find_ALE(graph, exonsDB, transcripts)
        if events:
            no_events[current_id] += 1
            write_GFF(events, exonsDB, no_events)

if __name__ == '__main__':
    main()
