"""
Usage:
    4c make-index [--flank=N] <genome> <site> <name>

Options:
    --flank=N  length of flag on left and right of
               each site to include in the index [default: 100]

Arguments:
    genome     fasta file of genome to process
    site       sequence of restriction site to use
    name       name of the restriction enzyme. Used as the
               name of the output directory for the index
               and as the prefix for the .info file.

Description:
    Create a bowtie index of sequences flanking selected
    restriction sites.
    
    Requires: bowtie and bowtie-build must be on PATH.
    Notes:
      1) bowtie2 is not yet supported.
      2) fragments with gaps >1kb are excluded from index
         but included in the .info file
      3) all output goes into directory <name>.  <name> must not
         exists.
      4) lmap and rmap are not yet implemented   

    For each restriction site found in the genome, sequences to the
    left and right are included in the index.  They are joined
    together with 4 Ns and named
   
    >enzyme_chrom_start0_end1
    
    In addition, a second file (<name>.info) contains details for each 
    fragment in the format
    
    enzyme_chrom_start0_end1|included|chrom|start0|end1|lmap|rmap
    
    where start0 is the 0-based start index (not including the site
    itself) and end1 is the 1-based end index such that end1 - start0
    is the actual length of the fragment.  Fragments that contain
    large gaps (>500nts) are not included in the index.  lmap and rmap
    are indicators [n|y|u] for the mappability of the left and right
    end, respectively.  They are strings of length [flank] such that
    lmap[i-1] represent the mappability of the +strand fragment of
    length i and rmap[i-1] represents the mappability of the right
    -strand fragment (i.e. starting from the right restriction
    site). n = not mappable; y=mappable; u = not tested
            
"""

import logging
import docopt
from .validators import is_valid_infile
import os
import sys
import numpy
import subprocess
from Bio import SeqIO

def validate(args):
    if not is_valid_infile(args["<genome>"]):
        logging.error("file '{<genome>}' not found".format(**args))
        sys.exit(1)
    args["<site>"] = args["<site>"].upper()
    if os.path.exists(args["<name>"]):
        logging.error("'{<name>}' already exists".format(**args))
        sys.exit(1)
    try:
        args["--flank"] = int(args["--flank"])
    except ValueError:
        logging.error("--flank requires an int value")
        sys.exit(1)
    
    logging.info("***** Creating new genome index *****")
    logging.info("Genome:     %s", args["<genome>"])
    logging.info("RE site:    %s", args["<site>"])
    logging.info("RE name:    %s", args["<name>"])
    logging.info("flank:      %d", args["--flank"])

    return args
    

def main(cmdline):
    args = validate(docopt.docopt(__doc__, argv=cmdline))
    make_index(args)


def make_index(genome, site, re_name, flank):
    out_dir = re_name
    os.mkdir(out_dir)
    
    out_fasta_fh = open(os.path.join(out_dir, "%s.fa" % re_name), "w")
    out_info_fh  = open(os.path.join(out_dir, "%s.info" % re_name), "w")
    make_fasta_file(genome, site, re_name, flank, 
            out_fasta_fh, out_info_fh)
    out_fasta_fh.close()
    out_info_fh.close()
    
    logging.info("Creating bowtie index")
    fnull = open(os.devnull, 'w')
    try:
        subprocess.check_call(["bowtie-build", 
            "%s.fa" % os.path.join(out_dir, args.name), 
            os.path.join(out_dir, args.name)],
            stdout = fnull, 
            stderr = fnull)
    except OSError:
        logging.exception("Could not find bowtie-build")
    except subprocess.CalledProcessError:
        logging.exception("bowtie-build did not finish correctly")
    logging.info("DONE")



def make_fasta_file(fasta, site, name, flank_len, out_fasta, out_info):
    fa_fmt   = ">{0}_{1}_{2}_{3}\n{4}NNNN{5}\n"
    info_fmt = "{0}_{1}_{2}_{3}|{4}|{1}|{2}|{3}|todo|todo\n"
    site     = site.upper()
    site_len = len(site)
    min_frag_len = flank_len + 2
    frag_lengths = []
    for rec in SeqIO.parse(fasta, "fasta"):
        pos      = 0
        seq      = rec.seq.upper()
        seq_find = seq.find
        seq_id   = rec.id
        logging.info("processing %s", seq_id)
        
        left_site = seq_find(site, start = pos)
        if left_site == -1:
            logging.warning("No site found in sequence %s", seq_id)
            break
        while True:
            right_site = seq_find(site, start = left_site + site_len)
            if right_site == -1:
                break
            start0   = left_site + site_len
            end1     = right_site
            frag_len = end1 - start0
            frag_lengths.append(min(frag_len, 10000))
            uncalled_n = seq.count("N", start = start0, end = end1)
            if uncalled_n < 1000 and frag_len > min_frag_len:
                out_fasta.write(fa_fmt.format(name, seq_id, start0, end1,
                    seq[start0:(start0 + flank_len)],
                    seq[(end1 - flank_len):end1]))
                out_info.write(info_fmt.format(name, seq_id, start0, end1,
                    "included"))
            elif uncalled_n >= 1000:
                out_info.write(info_fmt.format(name, seq_id, start0, end1,
                    "excluded[N=%d]" % uncalled_n))
                logging.debug("Fragment excluded due to Ns:      %s:%d-%d",
                        seq_id, start0, end1)
            else:
                out_info.write(info_fmt.format(name, seq_id, start0, end1,
                    "excluded[L=%d]" % frag_len))
                logging.debug("Fragment excluded short length:      %s:%d-%d (%d)",
                        seq_id, start0, end1, frag_len)

            left_site = right_site
    log_frag_len_histogram(frag_lengths, 500)


def log_frag_len_histogram(length_list, binsize):
    counts = numpy.bincount([x / binsize for x in length_list])
    bin_edges = ["[%5d - %5d[" % (i * binsize, (i + 1) * binsize) for i in range(len(counts))]
    bin_edges[-1] = ">%5d         " % (len(counts) * binsize)
    scale = 40.0 / max(counts) 
    logging.info("Size distribution of included fragments (mean = %d):",
            numpy.mean(length_list))
    for i in range(len(counts)):
        logging.info("%s |%s", bin_edges[i], int(counts[i] * scale) * "*")

