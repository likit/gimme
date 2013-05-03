===================================================
Gimme: A transcripts assembler based on alignments.
===================================================


Credits
-------

The program is developed in laboratory of genomics, evolution
and development (GED lab), Michigan State University.

:Web site: http://ged.msu.edu.

:Author: Likit Preeyanon, preeyano@msu.edu
:Advisor: C. Titus Brown, ctb@msu.edu

Copyright and license
---------------------

The prgram is Copyright Michigan State University.
The code is freely available for use and re-use under GNU GPL license.
See LICENSE.txt or http://www.gnu.org/licenses/.

Publication
-----------

Gimme is unpublished. A manuscript is in preparation.

Download
--------

Source code is available at https://github.com/ged-lab/gimme.git.

Installation
------------

Run python setup.py in the main directory to download and install required packages.

Running Gimme
-------------

Gimme should be able to run on any platform with Python 2.7 interpreter.
You can simply run::

    $ python ./src/gimme.py <input file>

Input
-----

Gimme can read an input file in PSL or BED format.
Use gff2bed.py in utils directory to convert GFF file to BED file.

Output
------

Output is written to standard output in BED format, which can be visualized
on UCSC genome browser or other browsers.

By default, gene models built by Gimme contain a minimum number of isoforms.
Use --max to force Gimme to report a minimum number of isoforms.
You can also use a script in utils to find a minimum set of transcripts.
See Utilities for more detail.

Example
-------

1. Assemble transcripts from sample data::

    $ python ./src/gimme.py sample_data/sample.psl > sample.bed

2. Obtain a maximum number of isoforms::

    $ python ./src/gimme.py -x sample_data/sample.psl > sample.max.bed

3. Run Gimme with multiple input files::

    $ python ./src/gimme.py sample1.psl sample2.psl sample3.psl > sample.all.bed

4. Run Gimme with user defined parameters::

    $ python ./src/gimme.py --min_utr=200 --max_intron=100000 --gap_size=15 sample.psl > sample.all.bed

5. See a program's help::

    $ python ./src/gimme.py -h or --help

Parameters
----------

GAP_SIZE, --gap_size=50
Introns smaller than GAP_SIZE) are filled to construct a more complete exon.

MAX_INTRON, --max_intron=300000
The maximum intron size (bp) allowed. A transcript is split into smaller parts
if it contains an intron longer than MAX_INTRON.

MIN_UTR, --min_utr=100
Alternative UTRs smaller than MIN_UTR are merged to overlapping exons.

MIN_TRANSCRIPT_LEN, --min_transcript_len=300
The minimum length (bp) for multiple exon transcript.

MIN_SINGLE_EXON_LEN, --min_single_exon_len=500
The minimum length (bp) for a single exon gene.

MAX_ISOFORMS, --max_isoforms=20
The maximum number of isoforms allowed without -x option.
Gimme searches for a minimum number of isoforms if the maximum number exceeds MAX_ISOFORMS.

-x, --max
Tell Gimme to search for report all putative isoforms.

--debug
Run Gimme with parameters set for debugging.

-v, --version
Print out a version number.

-h, --help
Print out a help message.

Running Tests
-------------

Run nosetests in the main directory to run all tests.

Utilities
---------

Gimme contains many useful utilities that work with PSL, BED and SAM format.
Some programs are useful for building gene models.
Others are useful for working with reads, assembly sequences etc.
