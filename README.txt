Gimme: a transcripts assembler based on alignments.
Author: Likit Preeyanon
Email: preeyano@msu.edu

Copyright and license
=====================
The prgram is Copyright Michigan State University.
The code is freely available for use and re-use under GNU GPL license.
See LICENSE.txt for more detail.

Publication
===========

Gimme is unpublished. A manuscript is in preparation.

Required package
====================
Gimme requires Networkx package to manipulate a graph data structure.
You can download Networkx from http://networkx.lanl.gov/download.html.
Please make sure networkx is in your PYTHONPATH or standard Python libraries
on your computer.

Download
========

Gimme is available on Github at git://github.com/likit/gimme.git.

Installation
============

No installation is required. However, you need to have Python 2.7 installed.

Running Gimme
=============

Gimme should be able to run on any platform with Python 2.7 interpreter.
You can simply run:

$ python gimme.py [input file].

Input
=====

Gimme reads input file in PSL format (BLAT and GMAP support this format).

Output
======

Output is written to standard output in BED format, which can be visualized
on UCSC genome browser or other browsers.

By default, gene models built by Gimme contain a maximum number of isoforms.
Use --min to force Gimme to report a minimum number of isoforms.
You can also use a script in utils to find a minimum set of transcripts.
See Utilities for more detail.

Example
=======

1.Assemble transcripts from sample data:

$ python ./src/gimme.py sample_data/sample.psl > sample.bed


2.Obtain a minimum number of isoforms:

$ python ./src/gimme.py --min sample_data/sample.psl > sample.bed

3.Run Gimme with multiple input files:

$ python ./src/gimme.py sample1.psl sample2.psl sample3.psl > sample.all.bed

4.Run Gimme with user defined parameters:

$ python ./src/gimme.py --MIN_UTR=200 --MAX_INTRON=100000 --GAP_SIZE=15 sample.psl > sample.all.bed

See a program's help or Parameters section for more detail.

5.See a program's help:

$ python ./src/gimme.py -h or --help

Parameters
==========

This version of Gimme fills up all small gaps in alignments.
The default gap is 21bp and is specified in GAP_SIZE parameter.

You can also modify a length of minimal UTR by specifying a MIN_UTR value.
The default is 100 bp. This value only is used as a cut off for alternative UTRs.

MAX_INTRON is a maximum size of an intron. Oversized introns will be excluded.

Running Tests
=============

Run nosetests in the main directory to run all tests in tests directory.

Utilities
=========

Gimme contains many useful utilities that work with PSL, BED and SAM format.
Some programs are useful for building gene models. Others are useful for
working with reads, assembly sequences etc.

get_min_path.py
---------------

The program finds a minimal set of transcripts that contains all edges from
an exon graph.

Example:

$ python ./src/utils/get_min_path.py models.bed > models.min.bed
