"""
Usage:
    4c fragcount <outdir> <bam> ...

Arguments:
    outdir  output directory
    bam     bam file(s) to process

Description:
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

    Note that that multiple processes are started in parallel and that
    therefore this should be run on a processing node.
"""




import docopt
from . import validators as val
import sys
import os
import logging
import collections
import subprocess
import shlex
import pysam

BUFSIZE = 81920

def check_args(args):
    schema = {"<outdir>": (lambda x: val.is_dir(x, allow_existing=True),
                           "{<outdir>} is not a valid target directory".format(**args)),
              "<bam>": (val.is_valid_bam_file_list,
                        "one or more of bam files were not found")}
    ok, errors = val.validate(args, schema)
    if not ok:
        for e in errors:
            logging.error(e)
        return None
    else:
        return args
    
def main(cmdline):
    args = docopt.docopt(__doc__, argv=cmdline)
    args = check_args(args)
    if args is not None:
        fragcount(args["<outdir>"],
                  args["<bam>"])
    else:
        sys.exit(1)
        
def fragcount(outdir, bam_lst):
    logging.info("output directory: %s", outdir)
    if not os.path.exists(outdir):
        logging.info("  creating...")
        os.mkdir(outdir, 0700)
    else:
        logging.info("  already exists")
    count(bam_lst, outdir)

def list2():
    return [0, 0]

class RunSummary(object):
    def __init__(self):
        self.frag_count  = {}
        self.valid_count = {}
    def count(self, is_ok, qname, rname, side):
        sample = qname.split(":")[0]
        if sample not in self.frag_count:
            self.frag_count[sample]  = collections.defaultdict(list2)
            self.valid_count[sample] = [0, 0]
        self.valid_count[sample][is_ok] += 1
        if is_ok == 1:
            self.frag_count[sample][rname][side] += 1
    def write(self, out_dir, all_frag_list):
        for sample, counts in self.valid_count.items():
            logging.info("%25s alignments: %8d valid, %8d invalid",
                    sample, counts[1], counts[0])
            row     = "{0}\t{1}\t{2}\t{3}\t{4}\t{5:.8f}\n"
        sorters = []
        for sample, counts in self.frag_count.items():
            out_fh = open(os.path.join(out_dir, sample + ".bed"), "w")
            sorter = subprocess.Popen(
                    shlex.split("sort -S1G -k1,1 -k2,2g"),
                    stdin     = subprocess.PIPE,
                    stdout    = out_fh,
                    bufsize   = BUFSIZE,
                    close_fds = True)
            total = self.valid_count[sample][1]
            for frag in counts:
                _, chrom, start0, end1 = frag.rsplit("_", 3)
                left_n, right_n = counts[frag]
                sorter.stdin.write(row.format(
                    chrom, start0, end1,
                    left_n, right_n,
                    float(left_n + right_n) / total))
            sorter.stdin.close()
            sorters.append(sorter)
       
        # create one bed file that contains all possible fragments
        out_fh = open(os.path.join(out_dir, "all_fragments.bed"), "w")
        sorter = subprocess.Popen(
                shlex.split("sort -S1G -k1,1 -k2,2g"),
                stdin     = subprocess.PIPE,
                stdout    = out_fh,
                bufsize   = BUFSIZE,
                close_fds = True)
        sorters.append(sorter)
        bed = "{0}\t{1}\t{2}\t{3}\n"
        for frag in all_frag_list:
            enz, chrom, start0, end1 = frag.rsplit("_", 3)
            sorter.stdin.write(bed.format(chrom, start0, end1, enz))
        sorter.stdin.close()
        logging.info("waiting for sort processes to finish")
        for sorter in sorters:    
            sorter.wait()
            if sorter.returncode != 0:
                logging.error("sort on output exited with error code %d",
                        sorter.returncode)
            else:
                logging.debug("sort on output exited normally")


#TODO:  if this is switched to a different aligner this function may have to be
#TODO:       modified
def count(bam_file_names, out_fh):
    frag = RunSummary()
    all_frag_list = None
    for bam_name in bam_file_names:
        bamf = pysam.Samfile(bam_name)
        if all_frag_list is None:
            all_frag_list = bamf.references
        ref_len = check_bam(bamf)
        for aln in bamf:
            if aln.is_unmapped:
                continue
            strand = "+"
            if aln.is_reverse:
                strand = "-"
            start0 = aln.pos
            end1   = aln.aend
            is_ok  = 0
            side   = None
            if strand == "+" and start0 == 0:
                is_ok = 1
                side  = 0
            elif strand == "-" and end1 == ref_len:
                is_ok = 1
                side  = 1
            frag.count(is_ok, aln.qname, bamf.getrname(aln.tid), side)
        bamf.close()
    frag.write(out_fh, all_frag_list)
    logging.info("DONE")

def check_bam(bamf):
    ref_lengths = list(set(bamf.lengths))
    if len(ref_lengths) == 1:
        logging.info("%s: ref sequence lengths: %d",
                bamf.filename, ref_lengths[0])
    else:
        logging.error("%s: ref sequences of different lengths found",
            bamf.filename)
        sys.exit(1)
    header = bamf.header.keys()
    if 'SQ' not in header:
        logging.error("%s: header is missing reference sequences (SQ)",
                bamf.filename)
        sys.exit(1)
    return ref_lengths[0]
