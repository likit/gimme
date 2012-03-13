'''The script reads gene models in BED format
and plot a histogram of transcript lengths

'''

import sys
import csv
import matplotlib.pyplot as plot

def get_len_dist(filename):
    print >> sys.stderr, "Reading from %s" % filename
    lengths = []
    reader = csv.reader(open(filename), dialect="excel-tab")
    for line in reader:
        lengths.append(sum([int(size) for size in line[10].split(",")]))

    return lengths

def main(filename):
    lengths = get_len_dist(filename)
    plot.hist(lengths, 200)
    plot.xlabel("length(bp)")
    plot.ylabel("transcript")
    plot.show()


if __name__=="__main__":
    main(sys.argv[1])
