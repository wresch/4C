Tools for 4C data analysis
==========================

As of now this project is not set up to work with
setuptools/distribute.  Instead, it is a self-contained directory
that can be directly executed with python, i.e.

    python 4c

Global help text can be displayed with

    python 4c -h

The tools are setup as subcommands

    python 4c {make-index,align,fragcount,bin}

and help for each can be obtained with

    python 4c subcommand -h

Subcommand descriptions (in brief):
-----------------------------------
Subcommands listed in order relevant for analyzing paired end solexa 4C
experiments:

<table>
<tr><td>make-index</td><td>Create bowtie index of restriction fragments (with flanking region)</td></tr>
<tr><td>align</td><td>Align 4c pairs to restriction fragment index</td></tr>
<tr><td>fragcount</td><td>Summarize the number of reads per detected fragment</td></tr>
<tr><td>bin</td><td>test help text</td></tr>
</table>

Dependencies
------------
###Python libraries:
* argparse (either python >= 2.7 or backported separate library)
* Biopython
* pybedtools
* pysam
* numpy

###External executables on path:
* bowtie and bowtie-build
* samtools
