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


def parse_reps(reps_file):
    kept = set()
    n = 0

    print >> sys.stderr, 'reading sequences from reps file'
    for line in open(reps_file):
        if line.startswith('>'):
            kept.add(line.strip().lstrip('>'))
            n += 1
            if n % 10000 == 0:
                print >> sys.stderr, '...', n
        else:
            continue

    return kept


def main(args):
    bed_file = args[1]
    reps_file = args[2]
    parseBED(bed_file, parse_reps(reps_file))


if __name__ == '__main__':
    main(sys.argv)
