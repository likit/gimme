'''The script counts reads mapped across each splice junction
from in gene models. It requires gene models in BED format
and read alignments in BAM format. Note that reads must be
mapped to transcripts derived from the same gene models.

Transcripts can be obtained using get_transcript_seq.py provided
in this package.

Note: count_spliced_reads2.py only counts reads that map to
splice junctions at the position away from either end of a read.

By default, the distance from either end of a read is DIST_END=5.

'''

import sys, csv
import pysam

DIST_END = 5

class Junction(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.coord = self.__str__()
        self.transcript_pos = None

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


def get_junction_position(row, junction):
    '''Returns a position of a junction on transcripts.

    '''

    junctions = []
    chrom = row[0]
    chrom_start = int(row[1])
    block_starts = [int(s) for s in row[-1].split(',')]
    block_sizes = [int(s) for s in row[-2].split(',')]

    for i in range(len(block_starts)):
        start = chrom_start + block_starts[i]
        end = block_sizes[i] + start

        try:
            next_start = chrom_start + block_starts[i+1]
        except IndexError:
            pass
        else:
            junc = '%s:%d-%d' % (chrom, end + 1, next_start - 1)
            junctions.append(junc)
    try:
        position = sum(block_sizes[:junctions.index(str(junction)) + 1])
    except ValueError:
        print >> sys.stderr, row[3], junction, 'not in list'
        return None
    else:
        return position


def get_junction(row):
    '''Retreive all junctions from a given transcripts.
    
    Argument:
        row : a list of all attributes of a transcript

    '''

    chrom = row[0]
    transcript = row[3]
    chrom_start = int(row[1])
    block_starts = [int(b) for b in row[-1].split(',')]
    block_sizes = [int(s) for s in row[-2].split(',')]

    for i in range(len(block_starts) - 1):
        start = block_starts[i] + block_sizes[i] + chrom_start + 1
        end = block_starts[i + 1] + chrom_start - 1
        junction = Junction(chrom, start, end)

        yield junction, transcript


def main(bedfile, samfile):
    '''Parses splice junctions from gene models in BED format.
    Then counts reads that mapped across a splice junction.
    Results are written to standard output.

    Arguments:
        bedfile : gene models in BED format
        samfile : alignments in BAM format

    '''

    junction_db = {}

    reader = csv.reader(open(bedfile), dialect='excel-tab')

    for n, row in enumerate(reader, start=1):
        for junction, transcript in get_junction(row):
            junction.transcript_pos = get_junction_position(row, junction)
            # print >> sys.stderr, junction.transcript_pos
            num_mapped_reads = 0
            for read in samfile.fetch(transcript,
                                                junction.transcript_pos,
                                                junction.transcript_pos + 1
                                                ):
                if not (junction.transcript_pos - read.pos < DIST_END or
                            read.aend - junction.transcript_pos < DIST_END):
                    num_mapped_reads += 1

            try:
                junction_db[str(junction)] += num_mapped_reads
            except KeyError:
                junction_db[str(junction)] = num_mapped_reads

        if (n % 1000) == 0:
            print >> sys.stderr, '...', n

    for junction, num_mapped_reads in junction_db.iteritems():
        print('{0}\t{1}'.format(junction, num_mapped_reads))


if __name__=='__main__':
    bedfile = sys.argv[1]
    samfile = pysam.Samfile(sys.argv[2], 'rb')
    main(bedfile, samfile)
