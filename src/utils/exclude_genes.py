'''The script excludes genes that have more than a given number of isoforms.

The default is 50.

'''

import sys
import csv


def parse_BED(filename):
    reader = csv.reader(open(filename), dialect='excel-tab')
    for row in reader:
        transcript_id = row[3].replace(':', '-')
        yield transcript_id, row


def main():
    infile = sys.argv[1]
    try:
        max_isoform = int(sys.argv[2])
    except IndexError:
        max_isoform = 50
        print >> sys.stderr, 'Default max isoform = 50'
    else:
        print >> sys.stderr, 'Max isoform = %d' % max_isoform

    excluded = 0
    current_id = None
    writer = csv.writer(sys.stdout, dialect='excel-tab')
    transcripts = []
    for transcript_id, row in parse_BED(infile):
        new_id = transcript_id.split('.')[0]
        if not current_id:  # first gene
            current_id = new_id
            transcripts.append(row)
        else:
            if new_id != current_id:
                if len(transcripts) <= max_isoform:
                    for trn in transcripts:
                        writer.writerow(trn)
                else:
                    excluded += 1

                current_id = new_id
                transcripts = []

            transcripts.append(row)

    if len(transcripts) <= max_isoform:
        for trn in transcripts:
            writer.writerow(trn)
    else:
        excluded += 1

    print >> sys.stderr, 'Total %d genes excluded.' % excluded

if __name__ == '__main__':
    main()
