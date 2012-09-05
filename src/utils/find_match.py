'''Reads BLAT output in Blast9 format and report

transcripts/genes matching sequences in a database.

(run blat with -out=blast9 to get an output in Blast9 format)

input :
    1) Gene models in BED format
    2) Alignments in Blast9 format

'''

import sys

''' a query must cover at least a given
percent of a sequence in a database.'''
PERCENT_COVERAGE = 80


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
        cols = line.split('\t')
        gene_id, trans_id = cols[3].strip().split('.')
        size = sum([int(block) for block in cols[10].split(',')])
        trans_len[trans_id] = size

    return trans_len

if __name__=='__main__':
    try:
        trans_len = get_transcript_length(sys.argv[1])
        genes, transcripts = find_match(sys.argv[2], trans_len)
    except IndexError:
        print >> sys.stderr, "find_match.py <BED file> <Blast9 file>"
    else:
        for gene in genes:
            print gene
