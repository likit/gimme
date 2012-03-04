'''The script reads an output file from
count_spliced_reads.py and plot cumulative
distribution of reads mapped to splice junctions.

count_spliced_reads.py is included in Gimme package.

'''
import sys

import matplotlib.pyplot as plot
from scipy.stats import cumfreq


def get_data(filename):
    '''Read read counts and return a dictionary

    containing the number of reads of each splice
    junction.

    '''
    junction_db = {}

    with open(filename, 'r') as infile:
        for line in infile:
            junction, reads = line.strip().split('\t')
            junction_db[junction] = int(reads)

    return junction_db


def main(filename):
    counts = get_data(filename)
    sorted_counts = sorted([v for v in counts.itervalues()])
    cumfreqs, lowlim, binsize, extrapoints = cumfreq(sorted_counts,
                                                        max(sorted_counts))
    norm_cumfreqs = cumfreqs / max(cumfreqs)
    plot.plot(norm_cumfreqs[:500], linewidth=1.5)
    plot.xlabel("mapped reads")
    plot.ylabel("splice junction")
    plot.show()


if __name__ == "__main__":
    main(sys.argv[1])
