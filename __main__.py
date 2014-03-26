"""

Usage:
    4c [-q|--quiet] <command> [<args>...]

Options:
    -q --quiet    Silence debug output

Commands (in order of use):
    make-index
    align
    fragcount
    bin
    mkfq
    
Description:
    A set of tools for generating raw per-fragment 4C counts
    from paired end reads where only the read spanning the primary
    restriction site is aligned to all possible restriction
    fragment ends from the genome. As of yet, no accounting for
    mappability is done (maybe in a later version).
    
Dependencies:
    Python
        - docopt
        - pybedtools
        - numpy
        - biopython
    Executables on path:
        - samtools
        - bowtie 1

Version: 0.1
"""

from __future__ import print_function
import docopt
import os
import sys
import logging

#TODO:  improve error handling in all modules (cleanup actions!!!)



BUFSIZE = 81920
################################################################################
# a outfile type for the argparser
################################################################################
def outfile_check(mode = "w", bufsize = BUFSIZE):
    """raises ArgumentTypeError if file s already exists; opens file for output
    otherwise"""
    def _outfile_check(s):
        if s == "-":
            return sys.stdout
        if os.path.exists(s):
            msg = "File %s already exists.  4c will not overwrite existing files" % s
            raise argparse.ArgumentTypeError(msg)
        return open(s, mode, bufsize)
    return _outfile_check
    
def infile_check(mode = "r", bufsize = BUFSIZE):
    """raises ArgumentTypeError if file s does not exist; handles '-' as special
    case for stdin"""
    def _infile_check(s):
        if s == "-":
            return sys.stdin
        if not os.path.exists(s):
            msg = "File %s does not exist" % s
            raise argparse.ArgumentTypeError(msg)
        if "w" in mode:
            raise ValueError("Infile mode can't contain 'w'")
        return open(s, mode, bufsize)
    return _infile_check

def outdir_check(s):
    if os.path.exists(s):
        raise argparse.ArgumentTypeError("Directory %s already exists", s)
    else:
        return s

def infile_name_check(s):
    if not os.path.exists(s):
        raise argparse.ArgumentTypeError("File %s does not exist", s)
    else:
        return s

################################################################################
# Root parser and setup
################################################################################

args = docopt.docopt(__doc__, options_first=True)
if args["--quiet"]:
    log_level = logging.INFO
else:
    log_level = logging.DEBUG
logging.basicConfig(
        level   = log_level,
        format  = "%(levelname)7s:%(funcName)s:%(lineno)d| %(message)s")

argv = [args['<command>']] + args['<args>']
if args["<command>"] == "make-index":
    from lib4c import make_index
    make_index.main(argv)
elif args["<command>"] == "align":
    from lib4c import align
    align.main(argv)
elif args["<command>"] == "fragcount":
    from lib4c import fragcount
    fragcount.main(argv)
elif args["<command>"] == "bin":
    from lib4c import bin
    bin.main(argv)
elif args["<command>"] == "mkfq":
    from lib4c import mkfq
    mkfq.main(argv)
else:
    print(__doc__, file=sys.stderr)
    sys.exit(1)
sys.exit()
