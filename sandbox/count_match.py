'''Count sequences that matched sequences in a database
from BLAST output in XML format.

'''

import sys

from Bio.Blast import NCBIXML

def find_match(input_file, evalue=0.01, verbose=True):
    '''Reads in Blast output in XML format and
    count all sequences that match sequences in a
    database with e-value less than a given threshold.

    '''

    handle = open(input_file)
    matches = set()
    blast_records = NCBIXML.parse(handle)
    for record in blast_records:
        query_name = record.query.split()[0]
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                subject = alignment.title.split()[0]
                if hsp.expect < evalue:
                    if verbose:
                        print '%s\t%e\t%s' % (query_name,
                                                hsp.expect,
                                                subject)
                    matches.add(query_name)
    return matches


def main(args):
    matches = find_match(args[1])
    print >> sys.stderr, 'Total match = %d' % len(matches)


if __name__=='__main__':
    main(sys.argv)
