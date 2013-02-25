'''The script reads BAM file and separate reads mapped with and without
mate into separate SAM files.

Reads are written in paired.bam and unpaired.bam accordingly.

Reads without a proper pair are ignored.

'''

import sys
import pysam


def split(infile):
    prefix = infile.rstrip(".bam")
    bamfile = pysam.Samfile(infile, 'rb')
    paired_file = pysam.Samfile(prefix + '_paired.sam',
                                        'w', template=bamfile)
    unpaired_file = pysam.Samfile(prefix + '_unpaired.sam',
                                        'w', template=bamfile)
    for read in bamfile.fetch():
        if read.is_unmapped:
            continue
        if read.is_paired:
            if read.mate_is_unmapped:
                unpaired_file.write(read)
            else:
                if read.is_proper_pair:
                    paired_file.write(read)

    paired_file.close()
    unpaired_file.close()


def main():
    infile = sys.argv[1]
    split(infile)


if __name__ == '__main__':
    main()
