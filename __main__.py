from __future__ import print_function
import argparse
import os
import sys
import logging
import textwrap

# pipeline modules
from lib4c.make_index import make_index
from lib4c.align      import align
from lib4c.fragcount  import fragcount
from lib4c.bin        import bin_frag
from lib4c.mkfq       import mkfq

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
# Root parser
################################################################################

cmdline  = argparse.ArgumentParser(
        description = " 4C PIPELINE ".center(70, "*"))
cmdline.add_argument("-q", "--quiet", action="store_true", default=False,
        help = "exclude debug messages")
commands = cmdline.add_subparsers(
        title       = "subcommands")

################################################################################
# interface for the make-index action
################################################################################
make_index_cmd = commands.add_parser("make-index",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        help            = "Create bowtie index of restriction fragments", 
        description     = textwrap.dedent("""\
            Create a bowtie index of sequences flanking selected restriction
            sites.
            
            Requires: bowtie and bowtie-build must be on PATH.
            Notes:    1) bowtie2 is not yet supported.
                      2) fragments with gaps >1kb are excluded from index
                         but included in the .info file
                      3) all output goes into directory <name>.  <name> must not
                         exists.
                      4) lmap and rmap are not yet implemented   

            For each restriction site found in the genome, sequences to the left 
            and right are included in the index.  They are joined together with 
            4 Ns and named  
           
            >enzyme_chrom_start0_end1
            
            In addition, a second file (<name>.info) contains details for each 
            fragment in the format
            
            enzyme_chrom_start0_end1|included|chrom|start0|end1|lmap|rmap
            
            where start0 is the 0-based start index (not including the site itself)
            and end1 is the 1-based end index such that end1 - start0 is the actual
            length of the fragment.  Fragments that contain large gaps (>500nts) are
            not included in the index.  lmap and rmap are indicators [n|y|u] for the
            mappability of the left and right end, respectively.  They are strings
            of length [flank] such that lmap[i-1] represent the mappability
            of the +strand fragment of length i and rmap[i-1] represents the
            mappability of the right -strand fragment (i.e. starting from the
            right restriction site). n = not mappable; y=mappable; u = not
            tested"""))
make_index_cmd.add_argument("genome", type=infile_check("r", BUFSIZE),
        help = "fasta file of genome to process")
make_index_cmd.add_argument("site",
        help = "restriction enzyme site sequence")
make_index_cmd.add_argument("name",
        help = "restriction enzyme name")
make_index_cmd.add_argument("--flank", type=int, default=100,
        help = """length of flag on left and right of each site to include in the
        index [%(default)s]""")
make_index_cmd.set_defaults(func = make_index)

################################################################################
# interface for the align action
################################################################################
align_cmd = commands.add_parser("align",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        help            = "Align 4c pairs to restriction fragment index", 
        description     = textwrap.dedent("""\
            Take a pair of compressed fastq files, determine which pairs are
            valid (i.e. have the 6-hitter and 4-hitter flank), extract the
            target sequence from the 6-hitter arm (not including the restriction
            site itself), and align the 6-hitter arm against a bowtie index
            containing only sequences flanking the 6-hitter sites.  Output is in
            bam format."""))
align_cmd.add_argument("read1", type = infile_name_check, 
        help = "First read of pair")
align_cmd.add_argument("read2", type = infile_name_check, 
        help = "Second read of pair")
align_cmd.add_argument("config", type = infile_check("r"), 
        help = "Configuration file describing flanks")
align_cmd.add_argument("genome",
        help = "bowtie indexed genome of restriction site flanks")
align_cmd.add_argument("-o", "--out", type = outfile_check("w", BUFSIZE),
        default="-", 
        help = "output file name [%(default)s]")
align_cmd.add_argument("-q", "--qual", 
        choices = ["solexa", "phred64", "phred33"],
        default = "phred33", 
        help    = "scale of quality values; passed on to bowtie [%(default)s]")
align_cmd.set_defaults(func = align)

################################################################################
# interface for the fragcount action
################################################################################
fragcount_cmd = commands.add_parser("fragcount",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        help            = "Summarize the number of reads per detected fragment",
        description     = textwrap.dedent("""\
            Take the bam output of the align action and generate a table of 
            read counts for each fragment per sample.  Multiple bam files
            can be run together if sample is spread over multiple alignment
            files.

            Output format:
                chrom|start0|end1|#left|#right|norm total
                | = tab

            where #left and #right are the number of reads mapping to the left
            and right sides of the fragment, respectively.  norm left/right are
            the normalized numbers (by aligned library size).

            The output directory will contain one file per sample plus one file
            called "all_fragments.bed" that contains a listing of all the
            restriction fragments parsed from the first bam header. This ***
            only *** makes sense if all the bam files were created with the same
            set of restriction fragments.
            """))
fragcount_cmd.add_argument("outdir", type=outdir_check,
        help = "output directory")
fragcount_cmd.add_argument("bam", nargs="+", type = infile_name_check, 
        help = "Bam file(s) to process")
fragcount_cmd.set_defaults(func = fragcount)

################################################################################
# interface for the bin action
################################################################################
#TODO: allow user to just specify size and step for bin generation based on 
#            genome coordinates.
bin_cmd = commands.add_parser("bin",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        help            = "test help text",
        description     = textwrap.dedent("""\
            Create bedgraph file for each fragcount output file showing the fraction of
            positive frgments per fixed-size bin.

            Currently user has to provide bins. 
            """))
bin_cmd.add_argument("outdir", type = outdir_check,
        help = "output directory")
bin_cmd.add_argument("bin_bed", type = infile_name_check, 
        help = "bed file listing bins used for analysis")
bin_cmd.add_argument("all_frag_bed", type = infile_name_check, 
        help = "bed file listing all fragments that were considered (created \
        by fragcount")
bin_cmd.add_argument("fragcount_bed", nargs="+", type = infile_name_check, 
        help = "fragcount file(s) to process")
bin_cmd.set_defaults(func = bin_frag)


################################################################################
# interface for the align action
################################################################################
mkfq_cmd = commands.add_parser("mkfq",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        help            = "Extract read pairs for samples (for example for submission to database)",
        description     = textwrap.dedent("""\
            Take a pair of compressed fastq files, determine which pairs are
            valid (i.e. have the 6-hitter and 4-hitter flank), and write the 
            pairs for each sample to separate files. Uses the same
            configuration file as the alignment command.
            """))
mkfq_cmd.add_argument("read1", type = infile_name_check, 
        help = "First read of pair")
mkfq_cmd.add_argument("read2", type = infile_name_check, 
        help = "Second read of pair")
mkfq_cmd.add_argument("config", type = infile_check("r"), 
        help = "Configuration file describing flanks")
mkfq_cmd.add_argument("-o", "--out", type = outdir_check,
        default = "split_fastq",
        help = "directory into which to save the split fastq files [%(default)s]")
mkfq_cmd.set_defaults(func = mkfq)



################################################################################
# run
################################################################################
args = cmdline.parse_args()
if args.quiet is True:
    log_level = logging.INFO
else:
    log_level = logging.DEBUG
logging.basicConfig(
        level   = log_level,
        format  = "%(levelname)-7s:%(asctime)s:%(funcName)s| %(message)s",
        datefmt = "%y%m%d:%H.%M.%S")
args.func(args)
