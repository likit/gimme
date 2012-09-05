'''The script reads alignments in PSL format
and cluster transcripts to loci.
Transcripts are written in separate files
based on loci.

'''

import sys
import csv

from collections import defaultdict, namedtuple
from bx.intervals.cluster import ClusterTree
from pygr import seqdb, sequtil

Alignment = namedtuple('Alignment', ['chrom', 'start', 'end', 'trans_id'])

def parse_psl(filename):
    reader = csv.reader(open(filename), dialect='excel-tab')
    
    for row in reader:
        chrom = row[13]
        trans_id = row[9]
        start = int(row[15])
        end = int(row[16])
        alignment = Alignment(chrom, start, end, trans_id)
        yield alignment


def build_cluster_trees(alignment_generator, cluster_distance=20):
    # arguments to ClusterTree are:
    # - Distance in basepairs for two reads to be in the same cluster;
    #   for instance 20 would group all reads with 20bp of each other
    # - Number of reads necessary for a group to be considered a cluster;
    #   2 returns all groups with 2 or more overlapping reads

    cluster_trees = defaultdict(lambda: ClusterTree(cluster_distance, 1))
    trans_db = {}
    trans_no = 0
    for alignment in alignment_generator:
        cluster_trees[alignment.chrom].insert(
                                                alignment.start,
                                                alignment.end,
                                                trans_no,
                                                #alignment.trans_id
                                                )
        trans_db[trans_no] = alignment
        trans_no += 1

    return dict(cluster_trees), trans_db


def get_sequence(cluster_trees, fasta_file, trans_db):
    print >> sys.stderr, 'reading %s' % fasta_file
    sequences = seqdb.SequenceFileDB(fasta_file, verbose=False)
    cluster_no = 0
    for chrom, cluster_tree in cluster_trees.items():
        for start, end, trans_nos in cluster_tree.getregions():
            print >> sys.stderr, '%d\t%s:%d-%d\t%d' \
                                % (cluster_no, chrom, start, end, len(trans_nos))
            with open('locus_%d.fa' % (cluster_no), 'w') as op:
                for no in trans_nos:
                    seqid = trans_db[no].trans_id
                    try:
                        sequence = sequences[seqid].seq
                    except KeyError:
                        pass
                    else:
                        sequtil.write_fasta(
                                            op,
                                            sequence,
                                            id=seqid,
                                        )
            cluster_no += 1


def main(argv):
    psl_file = argv[1]
    fasta_file = argv[2]

    cluster_trees, trans_db = build_cluster_trees(parse_psl(psl_file), 5)
    get_sequence(cluster_trees, fasta_file, trans_db)


if __name__=='__main__':
    main(sys.argv)
