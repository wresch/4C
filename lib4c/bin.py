"""
Usage:
    4c bin <outdir> <bins> <frags> <fragcount>...

Arguments:
    outdir      output directory
    bins        bins to use for binning in bed format
                (4 column bed format - bed3 + id)
    frags       all fragments that were considered in
                bed format (created by fragcount)
                also 4 column bed format.
    fragcount   fragment count file(s) (also created by
                fragcount - 6 columns)

Description:
    Create bed file for each fragcount output file with the following
    columns:
      - bin chrom 
      - bin start 
      - bin end 
      - total number of restriction fragments overlapping bin
      - sum of normalized 4C signal for bin
      - number of detected restriction fragments
      - fraction of detected restriction fragments (fracPos)

    Normalization is just by total readnumber. 

    Currently user has to provide bins. 


"""

import docopt
from . import validators as val
import os
import logging
import pybedtools as pbt
import sys

#TODO:  bad names for functions/variables - fix!

def check_args(args):
    schema = {"<outdir>": (lambda x: not os.path.exists(x),
                           "directory {<outdir>} already exists".format(**args)),
              "<bins>": (val.is_valid_infile,
                        "bin file {<bins>} does not exist".format(**args)),
              "<frags>": (val.is_valid_infile,
                         "frag info file {<frags>} does not exists".format(**args)),
              "<fragcount>": (val.is_valid_infile_list,
                              "one or more fragcount files were not valid")}
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
    if args is None:
        sys.exit(1)
    print(args)
    bin_frag(args["<outdir>"],
             args["<bins>"],
             args["<frags>"],
             args["<fragcount>"])

def bin_frag(outdir, bin_bed, all_frag_bed, fragcount_bed):
    logging.info("output directory: %s", outdir)
    logging.info("Bin file: %s", bin_bed)
    logging.info("All restriction fragments bed file: %s", all_frag_bed)
    logging.info("Number of fragdata files: %d", len(fragcount_bed))
    os.mkdir(outdir)

    # open bins file
    bins = pbt.BedTool(bin_bed)
    logging.info("read in %8d bins", len(bins))

    # open all frag file
    all_frag = pbt.BedTool(all_frag_bed)
    logging.info("read in %8d restriction fragments", len(all_frag))

    # match up bins with restriction fragments
    #TODO: stats on the result
    bins_with_any_frag = count_frags_per_bin(bins, all_frag)
    logging.info("bins that contained any fragments: %d", len(bins_with_any_frag))

    make_bedgraph_files(fragcount_bed, bins_with_any_frag, outdir)

    # cleanup
    pbt.cleanup(remove_all = True)


def count_frags_per_bin(bins, frags):
    """does not include empty bins in output"""
    bins_with_frag = bins.intersect(frags, wao = True)\
            .filter(intersect_by_overlap)\
            .groupby(g = [1, 2, 3], c = [9], o = ["count"])
    return bins_with_frag

def count_frags_per_bin_empty(bins, frags):
    """includes empty bins in output"""
    def _overlap_to_bool(x):
        if int(x[10]) > 0.51 * (int(x[6]) - int(x[5])):
            # overlap is > 50%
            x[10] = '1' # change overlap into 0/1
        else:
            x[10] = '0' # change overlap into 0/1
            x[9]  = '0' # don't include fragcount in this interval
        return x
    def _frac_pos_bg(x):
        fracPos = float(x[5]) / float(x[3])
        x.append(str(fracPos))
        return x
    # the arguments for groupby use 1-based indexing into the columns
    # very confusing
    bins_with_count = bins.intersect(frags, wao = True)\
            .each(_overlap_to_bool)\
            .groupby(g = [1, 2, 3, 4], c = [10,11], o = ["sum", "sum"])\
            .each(_frac_pos_bg)
    return bins_with_count


def intersect_by_overlap(x):
    o = int(x[8])
    if o == 0:
        return False
    frag_len = int(x[6]) - int(x[5])
    return o > 0.51 * frag_len


def make_bedgraph_files(samples, bins, outdir, color="166,86,40", smoothingWindow=5):
    for sample in samples:
        # name based on file name
        name = os.path.basename(sample).split(".")[0]
        outfile_name = os.path.join(outdir, name + ".bed")
        logging.info("Processing %s", name)

        # just count the number of positive fragments per bin
        count_frags_per_bin_empty(bins, pbt.BedTool(sample))\
                .saveas(outfile_name)
