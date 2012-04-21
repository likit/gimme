'''The script reads a sequence in FASTA format and
identify a repeat region due to misassembly (a bubble)
and collapse it.

Author: Likit Preeyanon
Email: preeyano@msu.edu
'''

import sys
from collections import namedtuple

import networkx as nx
import matplotlib.pyplot as plt

KMERSIZE = 27

Sequence = namedtuple('Sequence', ['id', 'seq'])

def parse_fasta(filename):
    seqid = None
    for line in open(filename):
        if line.startswith('>'):
            newid = line.strip().lstrip('>')
            if seqid == None:
                seqid = newid
                seq = ''
            elif newid != seqid:
                yield Sequence(seqid, seq)
                seqid = newid
                seq = ''
        else:
            seq += line.strip()

    yield Sequence(seqid, seq)


def construct_kmer(sequence, kmersize=KMERSIZE):
    if sequence.seq < kmersize:
        return None

    kmer_table = {}
    kmers = []
    path = []
    end = kmersize
    start = 0
    kmer_id = 0

    while len(sequence.seq) - start >= kmersize:
        kmer = sequence.seq[start:end]
        if kmer not in kmer_table:
            kmer_table[kmer] = kmer_id
            kmers.append(kmer)
            kmer_id += 1

        path.append(kmer_table[kmer])

        start += 1
        end += 1

    kmer_graph = nx.DiGraph()
    kmer_graph.add_path(path)
    print >> sys.stderr, 'total %d kmers' % len(kmers)

    return kmer_graph, kmers


def add_node(node_id, kmer_graph, path):
    succ = kmer_graph.successors(node_id)
    print node_id, succ
    if succ == []:
        return
    else:
        if len(succ) > 1:
            for s in succ:
                if s not in path:
                    path.append(s)
                    add_node(s, kmer_graph, path)
        else:
            for s in succ:
                path.append(s)
                add_node(s, kmer_graph, path)



def traverse(kmer_graph):
    '''Return a path from all nodes.'''
    root = kmer_graph.nodes()[0]
    path = [root]
    add_node(root, kmer_graph, path)

    return path


def collapse(kmer_graph, kmers):
    for bubble in nx.algorithms.cycles.simple_cycles(kmer_graph):
        mismatch = 0
        try:
            k2, k1 = bubble[-2], bubble[-1]
        except IndexError:
            continue
        else:
            for i in range(len(kmers[k1-1])):
                if kmers[k1-1][i] != kmers[k2][i]:
                    mismatch += 1

            if mismatch == 1:
                try:
                    kmer_graph.remove_node(k2)
                except nx.NetworkXError:
                    pass


def rebuild_sequence(path, kmers):
    sequence = kmers[0]
    for k in path[1:]:
        sequence += kmers[k][-1]

    return sequence

def main(argv):
    fasta_file = argv[1]

    try:
        kmersize = int(argv[2])
    except IndexError:
        kmersize = None

    #log_file = open(fasta_file+'.log', 'w')

    for n, sequence in enumerate(parse_fasta(fasta_file), start=1):
        print >> sys.stderr, sequence.id,
        kmer_graph, kmers = construct_kmer(sequence, kmersize)
        if nx.algorithms.cycles.simple_cycles(kmer_graph):
            collapse(kmer_graph, kmers)
            path = traverse(kmer_graph)
            new_sequence = rebuild_sequence(path, kmers)
            print '>%s\n%s' % (sequence.id, new_sequence)
        else:
            new_sequence = sequence.seq

        print >> sys.stderr, '%s\t%d\t%d'\
            % (sequence.id, len(sequence.seq) - len(new_sequence), len(sequence.seq))

        if n % 1000 == 0:
            print >> sys.stderr, '...', n

        #nx.draw(kmer_graph)
        #plt.show()

    #log_file.close()


if __name__=='__main__':
    main(sys.argv)
