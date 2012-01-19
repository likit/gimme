import sys
import csv
import pslparser

def parsePSL(filename):
    for pslObj in pslparser.read(open(filename)):
        strand = pslObj.attrib['strand']
        chrom = pslObj.attrib['tName']
        name = pslObj.attrib['qName']
        chromStart = pslObj.attrib['tStart']
        chromEnd = pslObj.attrib['tEnd']
        blockCount = pslObj.attrib['blockCount']
        blockSizes = ','.join([str(i) for i in pslObj.attrib['blockSizes']])
        blockStarts = ','.join([str(i - chromStart) for i in pslObj.attrib['tStarts']])
        thickStart = chromStart
        thickEnd = chromEnd

        writer = csv.writer(sys.stdout, dialect='excel-tab')
        writer.writerow((chrom,
                        chromStart,
                        chromEnd,
                        name,
                        255,
                        strand,
                        thickStart,
                        thickEnd,
                        '0,0,0',
                        blockCount,
                        blockSizes,
                        blockStarts))

if __name__=='__main__':
    parsePSL(sys.argv[1])
