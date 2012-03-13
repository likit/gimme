'''The script reads sequencing reads in fastq format
and split them into two separate files according to their
position in the file. The results are files containing
reads from odd and even positions respectively.

'''

import sys

def split(fastq_file, output_file1, output_file2):
    output_file1 = open(output_file1, 'w')
    output_file2 = open(output_file2, 'w')

    fastq_file = open(fastq_file)
    num_read = 0
    read_name = fastq_file.readline().strip()
    while read_name:
        read_seq = fastq_file.readline().strip()
        plus_sign = fastq_file.readline().strip()
        read_qual = fastq_file.readline().strip()

        output = output_file1 if (num_read % 2 == 0) else output_file2

        print >> output, '%s\n%s\n%s\n%s\n' % (
                                    read_name,
                                    read_seq,
                                    plus_sign,
                                    read_qual)
        num_read += 1
        if num_read % 10000 == 0:
            print >> sys.stderr, '...', num_read

        read_name = fastq_file.readline().strip()


if __name__=='__main__':
    split(sys.argv[1], sys.argv[2], sys.argv[3])
