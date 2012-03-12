'''The script reads alignements in PSL format
and filters out transcripts containing small
gaps (<MIN INTRON).

'''

import sys
import csv

MIN_INTRON = 70

reader = csv.reader(open(sys.argv[1]), dialect='excel-tab')
writer = csv.writer(sys.stdout, dialect='excel-tab')

removed = 0
for n, rows in enumerate(reader, start=1):
    chrom = rows[-8]
    transcript_id = rows[9]
    if transcript_id.find(chrom) < 0:
        removed += 1
    else:
        writer.writerow(rows)

    if n % 10000 == 0:
        print >> sys.stderr, '...', n, 'removed', removed
