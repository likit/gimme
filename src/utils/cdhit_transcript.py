'''The script reads output from CDHIT-EST and selects
a representative sequence of each cluster and write
it to standard output in BED format.

'''

import sys
import csv

def parseBED(bed_file, kept):
    reader = csv.reader(open(bed_file), dialect='excel-tab')
    writer = csv.writer(sys.stdout, dialect='excel-tab')

    print >> sys.stderr, 'filtering sequences...'
    for n, line in enumerate(reader, start=1):
        tr_name = line[3]
        if tr_name in kept:
            writer.writerow(line)

        if n % 10000 == 0:
            print >> sys.stderr, '...', n


def parse_clstr(clstr_file):
    kept = set()

    print >> sys.stderr, 'reading sequences from .clstr file'
    for n, line in enumerate(open(clstr_file), start=1):
        if '*' in line:
            tr_name = line.split()[2]
            tr_name = tr_name.lstrip('>').rstrip('...')
            kept.add(tr_name)

        if n % 10000 == 0:
            print >> sys.stderr, '...', n

    return kept


def main(args):
    bed_file = args[1]
    clstr_file = args[2]
    parseBED(bed_file, parse_clstr(clstr_file))


if __name__=='__main__':
    main(sys.argv)
