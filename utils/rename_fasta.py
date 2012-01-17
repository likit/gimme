#!/opt/local/bin/python
''' This script rename sequence names in fasta format.
The new name will be in a format of
>directory name + a transcript number.
For example, >assembly1-39.
Output is written to stdout.

Author  Likit Preeyanon.
Email   preeyano@msu.edu

'''

import sys, os

input_file = sys.argv[1]
n = 0
dirname = os.getcwd().split(os.path.sep)[-1]

for line in open(input_file):
    if line.startswith('>'):
        print >> sys.stdout, '>%s-%d' % (dirname, n)
        n += 1
    else:
        print >> sys.stdout, line.strip(os.linesep)
