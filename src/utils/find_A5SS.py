'''The script identifies alternative 5' splice site from a BED file.
Output is written in GFF format suitable for differential
exon usage analysis using MISO.

'''

import sys
import csv

import networkx as nx
from bx.intervals.intersection import IntervalTree


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


class Intron(object):
    def __init__(self, chrom, start, end, geneID, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.geneID = geneID
        self.strand = strand
        self.flanking_exons = set()

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


def get_introns(exons, intronsDB):
    introns = []
    for i in range(len(exons) - 1):
        ex1, ex2 = exons[i], exons[i + 1]
        intron_start = ex1.end + 1
        intron_end = ex2.start - 1
        intron = Intron(ex1.chrom, intron_start, intron_end,
                        ex1.geneID, ex1.strand)

        if str(intron) not in intronsDB:
            intronsDB[str(intron)] = intron

        intron.flanking_exons.add((ex1, ex2))
        introns.append(intron)
    return introns


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


def find_A5SS(graph, interval, exonsDB):
    exons = graph.nodes()
    exons = sorted([exonsDB[exon] for exon in exons], key=lambda x: x.start)
    altss_events = []
    altss_exons = []
    curr_exon = exons[0]
    down_exons = set(graph.successors(str(curr_exon)))
    up_exons = set([curr_exon])
    i = 0
    while i < (len(exons) - 1):
        next_exon = exons[i + 1]
        if (curr_exon.start == next_exon.start and
                curr_exon.end != next_exon.end):
            down_exons.intersection_update(graph.successors(str(next_exon)))
            up_exons.add(next_exon)
        else:
            if (len(down_exons) > 0 and len(up_exons) > 1):
                for up in up_exons:
                    for dn in down_exons:
                        altss_exons.append((up, exonsDB[dn]))
                        # print >> sys.stderr, "%s, %s" % (str(up), str(dn))

                # print >> sys.stderr, '*' * 40
                altss_events.append(altss_exons)

            curr_exon = next_exon
            down_exons = set(graph.successors(str(curr_exon)))
            up_exons = set([curr_exon])
            altss_exons = []
        i += 1
    return altss_events


def write_GFF(events, exonsDB, no_events, redundant):
    all_exons = set()
    for event in events:
        for exon in event:
            all_exons.add(exon)
    all_exons = sorted(list(all_exons), key=lambda x: x.end)

    first_exon = all_exons[0]
    last_exon = all_exons[-1]

    mrnaid = 1
    event_no = str(no_events[first_exon.geneID])
    geneID = first_exon.geneID + '.ev' + event_no
    output_list = []
    output = "%s\tA5SS\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s" % (
            first_exon.chrom, first_exon.start, last_exon.end,
            first_exon.strand, geneID, first_exon.geneID)
    output_list.append(output)

    unique_event = 0
    for event in events:
        event_exons = sorted(event, key=lambda x: x.end)
        first_exon = event_exons[0]
        last_exon = event_exons[-1]
        output = "%s\tA5SS\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.%d;Parent=%s" % (
                        first_exon.chrom, first_exon.start, last_exon.end,
                        first_exon.strand, geneID, mrnaid, geneID)
        output_list.append(output)

        exonid = 1
        first_exon = event_exons[0]
        last_exon = event_exons[-1]
        altss = "%s-%s" % (str(first_exon), str(last_exon))
        for exon in event_exons:
            output = "%s\tA5SS\texon\t%d\t%d\t.\t%s" + \
                            "\t.\tID=%s.%d.%d;Parent=%s.%d" \
                            % (exon.chrom, exon.start, exon.end,
                                exon.strand, geneID, mrnaid, exonid,
                                geneID, mrnaid)
            exonid += 1
            output_list.append(output)

        if altss not in redundant:
            redundant.add(altss)
            unique_event += 1
            mrnaid += 1

    if unique_event:
        for output in output_list:
            print output


def main():
    redundant = set()
    no_events = {}
    exonsDB = {}
    intronsDB = {}
    infile = sys.argv[1]
    graph = nx.DiGraph()
    current_id = None
    intron_interval = IntervalTree()
    for exons, transcript_id in get_exon_node(infile):
        new_id = transcript_id.split('.')[0]
        if not current_id:  # first gene
            for e in exons:
                exonsDB[str(e)] = e
            graph.add_path([str(e) for e in exons])
            introns = get_introns(exons, intronsDB)

            assert len(exons) == len(introns) + 1

            for intron in introns:
                intron_interval.insert_interval(intron)

            current_id = new_id
            no_events[current_id] = 0
        else:
            if new_id != current_id:
                if len(graph.nodes()) > 1:
                    for events in find_A5SS(graph, intron_interval, exonsDB):
                        no_events[current_id] += 1
                        write_GFF(events, exonsDB, no_events, redundant)

                graph = nx.DiGraph()
                exonsDB = {}
                intronsDB = {}
                current_id = new_id
                no_events[current_id] = 0
                intron_interval = IntervalTree()

            for e in exons:
                exonsDB[str(e)] = e
            graph.add_path([str(e) for e in exons])
            introns = get_introns(exons, intronsDB)

            assert len(exons) == len(introns) + 1

            for intron in introns:
                intron_interval.insert_interval(intron)

    if len(graph.nodes()) > 1:
        for events in find_A5SS(graph, intron_interval, exonsDB):
            no_events[current_id] += 1
            write_GFF(events, exonsDB, no_events, redundant)

if __name__ == '__main__':
    main()
