#!/opt/local/bin/python
''' This script rename sequence names in fasta format.
The new name will be in a format of
>prefix + a transcript number.
For example, >assembly1-39.
Output is written to stdout.

Author  Likit Preeyanon.
Email   preeyano@msu.edu

'''

import sys
import os

input_file = sys.argv[1]
try:
    prefix = sys.argv[2]
except:
    prefix = ''
else:
    prefix = prefix + '-'

n = 0

for line in open(input_file):
    if line.startswith('>'):
        print >> sys.stdout, '>%s%d' % (prefix, n)
        n += 1
    else:
        print >> sys.stdout, line.strip(os.linesep)
