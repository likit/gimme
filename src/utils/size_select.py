'''Selects a fasta sequence with longer than a given size.'''

import sys
from Bio import SeqIO

def select(inputFile):
    for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
        if len(seq_record) > MIN_SIZE:
            yield seq_record

if __name__=='__main__':
    try:
        inputFile = sys.argv[1]
        MIN_SIZE = int(sys.argv[2])
    except:
        print >> sys.stderr, "Usage: python sizeSelect.py <inFile> <size>"
    else:
        SeqIO.write(select(inputFile), sys.stdout, "fasta")
