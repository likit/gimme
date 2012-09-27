#! /usr/local/bin/python
'''The script parses PSL file (i.e. from BLAT output).
read method returns each alignment stored in an PSL object.

'''

import sys

class PSL(object):
    def __init__(self, **kwargs):
        self.attrib = kwargs

    def __getattr__(self, key):
        return self.attrib[key]


def read(fobj, comment=None):
    '''Return an object of an alignment.

    fobj = file object.

    comment = read() will ignore the line starting
    with comment character.

    '''
    n = 0
    for line in fobj:
        if comment and line.startswith(comment):
                continue

        attrib = {}
        rows = line.split()

        try:
            assert len(rows) == 21
        except AssertionError:
            print >> sys.stderr, '>%d' % n, line
            n += 1
            continue

        attrib['matches'] = rows[0]
        attrib['misMatches'] = rows[1]
        attrib['repMatches'] = rows[2]
        attrib['nCount'] = int(rows[3])
        attrib['qNumInsert'] = rows[4]
        attrib['qBaseInsert'] = rows[5]
        attrib['tNumInsert'] = rows[6]
        attrib['tBaseInsert'] = rows[7]
        attrib['strand'] = rows[8]
        attrib['qName'] = rows[9]
        attrib['qSize'] = int(rows[10])
        attrib['qStart'] = int(rows[11])
        attrib['qEnd'] = int(rows[12])
        attrib['tName'] = rows[13]
        attrib['tSize'] = int(rows[14])
        attrib['tStart'] = int(rows[15])
        attrib['tEnd'] = int(rows[16])
        attrib['blockCount'] = int(rows[17])
        attrib['blockSizes'] = [int(i) for i in rows[18].split(',')[:-1]]
        attrib['qStarts'] = [int(i) for i in rows[19].split(',')[:-1]]
        attrib['tStarts'] = [int(i) for i in rows[20].split(',')[:-1]]
        
        pobj = PSL(**attrib) # pobj = PSL object
        if attrib['blockCount'] == 0:
            print >> sys.stderr, '>%d' % n, line
            n += 1
            continue
        else:
            n += 1
            yield pobj

if __name__ == '__main__':
    pass
