import sys
import csv

MAX_ISOFORM = 50

def parse_BED(filename):
    reader = csv.reader(open(filename), dialect='excel-tab')
    for row in reader:
        transcript_id = row[3].replace(':', '-')
        yield transcript_id, row

def main():
    infile = sys.argv[1]
    current_id = None
    writer = csv.writer(sys.stdout, dialect='excel-tab')
    transcripts = []
    for transcript_id, row in parse_BED(infile):
        new_id = transcript_id.split('.')[0]
        print >> sys.stderr, current_id, new_id
        if not current_id: # first gene
            current_id = new_id
            transcripts.append(row)
        else:
            if new_id != current_id:
                if len(transcripts) <= MAX_ISOFORM:
                    for trn in transcripts: writer.writerow(trn)

                current_id = new_id
                transcripts = []

            transcripts.append(row)

    if len(transcripts) <= MAX_ISOFORM:
        for trn in transcripts: writer.writerow(trn)

if __name__=='__main__':
    main()
