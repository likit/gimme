'''Split merged paired reads from diginormed reads into two files.

Reads without a proper mate are written in single.fasta.

'''

import sys
import gzip
import os
from collections import namedtuple

Read = namedtuple('Read', ['name', 'seq'])


def split(input_file):

    if input_file.endswith('gz') or input_file.endswith('gzip'):
        fp = gzip.open(input_file)
    else:
        fp = open(input_file)

    output_filename = os.path.basename(input_file)

    op1 = open(output_filename + '_1', 'w')
    op2 = open(output_filename + '_2', 'w')
    single = open(output_filename + '_unpaired', 'w')

    F = None  # a forward read
    while True:
        try:
            if not F:
                '''For the first time, read F and R from the file.'''

                F = Read(fp.readline().strip(),
                        fp.readline().strip())
                R = Read(fp.readline().strip(),
                        fp.readline().strip())
            else:
                R = Read(fp.readline().strip(),
                        fp.readline().strip())
            if R.name == '':
                print >> single, "%s\n%s" % (F.name, F.seq)
                break
        except:
            raise
            break
        else:
            if F.name.rsplit("/")[0] == R.name.rsplit("/")[0]:
                if F.name.endswith("/1") and R.name.endswith("/2"):
                    print >> op1, "%s\n%s" % (F.name, F.seq)
                    print >> op2, "%s\n%s" % (R.name, R.seq)
                else:
                    # Swapped mates
                    print >> op1, "%s\n%s" % (R.name, R.seq)
                    print >> op2, "%s\n%s" % (F.name, F.seq)
            else:
                print >> single, "%s\n%s" % (F.name, F.seq)
                F = R

    op1.close()
    op2.close()
    single.close()

if __name__ == '__main__':
    input_file = sys.argv[1]
    split(input_file)
