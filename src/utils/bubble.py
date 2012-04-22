'''The script reads a sequence in FASTA format and
identify a repeat region due to misassembly (a bubble)
and collapse it.

Author: Likit Preeyanon
Email: preeyano@msu.edu
'''

import sys
from collections import namedtuple

import networkx as nx

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

    kmer_table = {} # stores all kmers with their node IDs.
    kmers = [] # stores all kmers in order as in a graph.
    path = [] # used to build a directed kmer graph

    end = kmersize
    start = 0

    '''Kmer id is unique for each kmer.
    A path contains a kmer id, not a node ID so
    that one can identify a loop in a graph.

    '''
    kmer_id = 0

    '''Node id is not the same as kmer id.
    A kmer can have one or more node id.
    It is used to indicate where a location of
    a kmer a graph.

    '''
    node_id = 0

    while len(sequence.seq) - start >= kmersize:
        kmer = sequence.seq[start:end]
        if kmer not in kmer_table:
            kmer_table[kmer] = [node_id]
            kmer_id += 1
        else:
            kmer_table[kmer].append(node_id)

        path.append(kmer_table[kmer][0])
        kmers.append(kmer)

        start += 1
        end += 1
        node_id += 1

    kmer_graph = nx.DiGraph()
    kmer_graph.add_path(path)
    print >> sys.stderr, 'total %d kmers' % len(kmers)

    return kmer_graph, kmers, kmer_table


def add_node(node_id, kmer_graph, path):
    succ = kmer_graph.successors(node_id)
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


def collapse(kmer_graph, kmers, kmer_table):
    discard = False
    removed = []
    for node in kmer_graph.nodes():
        preds = kmer_graph.predecessors(node)
        succ = kmer_graph.successors(node)
        if (len(succ) == 2 and discard):
            discard = False
        elif discard:
            kmer = kmers[node]
            removed.append(kmer_table[kmer][-1])
        else:
            if len(preds) == 2:
                mismatch = 0
                pre1 = kmers[preds[0]]
                pre2 = kmers[preds[1]]
                for i in range(len(pre1)):
                    if pre1[i] != pre2[i]:
                        mismatch += 1

                if mismatch == 1:
                    discard = True
    return removed


def rebuild_sequence(removed, kmers):
    keep = []
    for i in range(len(kmers)):
        if i not in removed:
            keep.append(kmers[i])

    sequence = keep[0]
    for i in range(len(keep[1:])):
        sequence += keep[i][-1]

    return sequence

def main(argv):
    fasta_file = argv[1]

    try:
        kmersize = int(argv[2])
    except IndexError:
        kmersize = None

    log_file = open(fasta_file+'.log', 'w')

    for n, sequence in enumerate(parse_fasta(fasta_file), start=1):
        print >> sys.stderr, sequence.id,
        kmer_graph, kmers, kmer_table = construct_kmer(sequence, kmersize)
        '''
        for node in kmer_graph.nodes():
            if len(kmer_graph.predecessors(node)) > 1:
                print node,\
                        kmers[node],\
                        kmer_table[kmers[node]],\
                        kmer_graph.predecessors(node)
        '''
        removed = collapse(kmer_graph, kmers, kmer_table)
        new_sequence = rebuild_sequence(removed, kmers)
        if removed:
            #print '>%s\n%s' % (sequence.id, sequence.seq)
            print '>%s\n%s' % (sequence.id+'_new', new_sequence)

        print >> log_file, '%s\t%d\t%d'\
            % (sequence.id, len(sequence.seq) - len(new_sequence), len(sequence.seq))

        #path = traverse(kmer_graph)
        #new_sequence = rebuild_sequence(path, kmers)
        #print '>%s\n%s' % (sequence.id, new_sequence)
        #new_sequence = sequence.seq

        if n % 1000 == 0:
            print >> sys.stderr, '...', n

        #nx.draw(kmer_graph)
        #plt.show()

    log_file.close()


if __name__=='__main__':
    main(sys.argv)
