'''Reads BLAT output in Blast9 format and report

transcripts/genes matching sequences in a database.

(run blat with -out=blast9 to get an output in Blast9
format)

input :
    1) transcript lengths (from get_transcript_length.py)
    2) alignments in Blast9 format

'''

import sys

''' a query must cover at least a given
percent of a sequence in a database.'''
PERCENT_COVERAGE = 70


def find_match(filename, trans_len):
    genes = set([])
    transcripts = set([])

    for line in open(filename):
        if line.startswith('#'):
            continue
        else:
            cols = line.strip().split('\t')
            gene, trans = cols[1].split('.')
            align_len = int(cols[3])
            if ((align_len/float(trans_len[trans]) * 100.0)
                            > PERCENT_COVERAGE):
                    genes.add(gene)
                    transcripts.add(trans)
    return genes, transcripts


def get_transcript_length(filename):
    trans_len = {}
    for line in open(filename):
        gene, trans, size = line.strip().split('\t')
        trans_len[trans] = int(size)

    return trans_len

if __name__=='__main__':
    trans_len = get_transcript_length(sys.argv[1])
    genes, transcripts = find_match(sys.argv[2], trans_len)
