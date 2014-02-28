"""
Usage:
    4c bin <outdir> <bins> <frags> <fragcount>

Arguments:
    outdir      output directory
    bins        bins to use for binning in bed format
    frags       all fragments that were considered in
                bed format (created by fragcount)
    fragcount   fragment count file (also created by
                fragcount

Description:
    Create bedgraph file for each fragcount output file showing the fraction of
    positive frgments per fixed-size bin.

    Currently user has to provide bins. 


"""

import docopt
import os
import logging
import pybedtools as pbt

#TODO:  bad names for functions/variables - fix!

def main(cmdline):
    args = docopt.docopt(__doc__, argv=cmdline)
    print(args)

def bin_frag(args):
    logging.info("output directory: %s", args.outdir)
    logging.info("Bin file: %s", args.bin_bed)
    logging.info("All restriction fragments bed file: %s", args.all_frag_bed)
    logging.info("Number of fragdata files: %d", len(args.fragcount_bed))
    data_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(data_dir, "data")
    logging.info("data dir: %s", data_dir)

    os.mkdir(args.outdir)

    # open bins file
    bins = pbt.BedTool(args.bin_bed)
    logging.info("read in %8d bins", len(bins))

    # open all frag file
    all_frag = pbt.BedTool(args.all_frag_bed)
    logging.info("read in %8d restriction fragments", len(all_frag))

    # match up bins with restriction fragments
    #TODO: stats on the result
    bins_with_any_frag = count_frags_per_bin(bins, all_frag)
    logging.info("bins that contained any fragments: %d", len(bins_with_any_frag))

    make_bedgraph_files(args.fragcount_bed, bins_with_any_frag, args.outdir)

    # cleanup
    pbt.cleanup(remove_all = True)


def count_frags_per_bin(bins, frags):
    """does not include empty bins in output"""
    bins_with_frag = bins.intersect(frags, wao = True)\
            .filter(intersect_by_overlap)\
            .groupby(grp = [1, 2, 3], opCols = [9], ops = ["count"])
    return bins_with_frag

def count_frags_per_bin_empty(bins, frags):
    """includes empty bins in output"""
    def _filter(x):
        o = int(x[10])
        if o == 0:
            return True
        elif o > 0.51 * (int(x[6]) - int(x[5])):
            return True
        else:
            x[10] = '0'
            return True
    def _overlap_to_bool(x):
        o = int(x[10])
        if o > 0:
            x[10] = '1'
        else:
            # fix normcount column for bins without positive frags
            x[9] = '0'
        return x
    def _frac_pos_bg(x):
        fracPos = float(x[5]) / float(x[3])
        x.name = str(fracPos)
        return x
    bins_with_count = bins.intersect(frags, wao = True)\
            .filter(_filter)\
            .each(_overlap_to_bool)\
            .groupby(grp = [1, 2, 3, 4], opCols = [10,11], ops = ["sum", "sum"])\
            .each(_frac_pos_bg).cut([0, 1, 2, 3])
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
        outfile_name = os.path.join(outdir, name + ".bg")
        logging.info("Processing %s", name)
        trackline = " ".join(
            ["track type=bedGraph",
             "name='%s' color=%s" % (name, color), 
             "visibility=full alwaysZero=on",
             "smoothingWindow=%d maxHeightPixels=60:80:120" % smoothingWindow])

        # just count the number of positive fragments per bin
        count_frags_per_bin_empty(bins,
                pbt.BedTool(sample)).saveas(outfile_name, trackline = trackline)
        foo = count_frags_per_bin_empty(bins,
                pbt.BedTool(sample))
