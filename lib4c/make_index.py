"""
Usage:
    4c make-index [--flank=N] [--site2=SEQ] <genome> <site> <name>

Options:
    --flank=N    length of flag on left and right of
                 each site to include in the index [default: 100]
    --site2=SEQ  if given, records for each fragment the distance
                 between the primary sites on the left and the right
                 and the secondary site targeted in the second
                 digestion. If no secondary site is present, the
                 numbers given are -1,-1 (these fragments are
                 termed "blind"). [default: ]
    
Arguments:
    genome       fasta file of genome to process
    site         sequence of restriction site to use
    name         name of the restriction enzyme. Used as the
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
    
    enzyme_chrom_start0_end1|included|chrom|start0|end1|lmap|rmap|site2_dist
    
    where start0 is the 0-based start index (not including the site
    itself) and end1 is the 1-based end index such that end1 - start0
    is the actual length of the fragment.  Fragments that contain
    large gaps (>20% of fragment length) are not included in the
    index.  lmap and rmap are indicators [n|y|u] for the mappability
    of the left and right end, respectively.  They are strings of
    length [flank] such that lmap[i-1] represent the mappability of
    the +strand fragment of length i and rmap[i-1] represents the
    mappability of the right -strand fragment (i.e. starting from the
    right restriction site). n = not mappable; y=mappable; u = not
    tested [this is not implemented yet].
            
"""

import logging
import docopt
from . import validators as val
import os
import sys
import numpy
import subprocess
from Bio import SeqIO

def check_args(args):
    schema = {"<genome>": (val.is_valid_fa_file,
                           "file {} not found or not a fasta file".format(args["<genome>"])),
              "<name>": (lambda x: not os.path.exists(x),
                         "'{}' already exists".format(args["<name>"])),
              "<site>": (val.is_seq, "'{}' is not a valid re site".format(args["<site>"])),
              "--flank": (val.is_number, "flank needs to be a number"),
              "--site2": (lambda x: x=="" or val.is_seq(x),
                          "'{}' is not a valid re site".format(args["--site2"]))}
    ok, errors = val.validate(args, schema)
    if not ok:
        for e in errors:
            logging.error(e)
        return None
    args["<site>"] = args["<site>"].upper()
    args["--site2"] = args["--site2"].upper()
    args["--flank"] = int(args["--flank"])
    return args

def main(cmdline):
    args = check_args(docopt.docopt(__doc__, argv=cmdline))
    if args is None:
        sys.exit(1)
    logging.info("***** Creating new genome index *****")
    logging.info("Genome:     %s", args["<genome>"])
    logging.info("RE site:    %s", args["<site>"])
    logging.info("secondary RE site: %s", args["--site2"])     
    logging.info("RE name:    %s", args["<name>"])
    logging.info("flank:      %d", args["--flank"])
    try:
        bowtie_version=subprocess.check_output(["bowtie", "--version"]).split("\n")[0]
        logging.info("bowtie:   %s", bowtie_version)
    except OSError:
        logging.error("bowtie executable not found on path")
        sys.exit(1)

    make_index(args["<genome>"],
               args["<site>"],
               args["<name>"],
               args["--flank"],
               args["--site2"])

def make_index(genome, site, re_name, flank, site2):
    out_dir = re_name
    os.mkdir(out_dir)
    
    out_fasta_fh = open(os.path.join(out_dir, "%s.fa" % re_name), "w")
    out_info_fh  = open(os.path.join(out_dir, "%s.info" % re_name), "w")
    make_fasta_file(genome, site, re_name, flank, 
            out_fasta_fh, out_info_fh, site2)
    out_fasta_fh.close()
    out_info_fh.close()
    
    logging.info("Creating bowtie index")
    fnull = open(os.devnull, 'w')
    try:
        subprocess.check_call(["bowtie-build", 
            "%s.fa" % os.path.join(out_dir, re_name), 
            os.path.join(out_dir, re_name)],
            stdout = fnull, 
            stderr = fnull)
    except subprocess.CalledProcessError:
        logging.exception("bowtie-build did not finish correctly")
    logging.info("DONE")

def make_fasta_file(fasta, site, name, flank_len, out_fasta, out_info, site2):
    fa_fmt    = ">{0}_{1}_{2}_{3}\n{4}NNNN{5}\n"
    info_fmt  = "{0}_{1}_{2}_{3}|{4}|{1}|{2}|{3}|todo|todo|{5}\n"
    site      = site.upper()
    site_len  = len(site)
    site2     = site2.upper()
    site2_len = len(site2)
    min_frag_len = flank_len + 2
    frag_lengths = []
    max_n_freq = 0.2
    for rec in SeqIO.parse(fasta, "fasta"):
        pos      = 0
        seq      = rec.seq.upper()
        seq_find = seq.find
        seq_id   = rec.id
        logging.info("processing %s", seq_id)
        
        left_site = seq_find(site, start = pos)
        if left_site == -1:
            logging.error("No site found in sequence %s", seq_id)
            sys.exit(1)
        while True:
            right_site = seq_find(site, start = left_site + site_len)
            if right_site == -1:
                break
            start0   = left_site + site_len #not including the RE site
            end1     = right_site #not including the RE site
            frag_len = end1 - start0
            frag_lengths.append(min(frag_len, 10000))
            uncalled_n = float(seq.count("N", start = start0, end = end1)) / max(frag_len, 1)
            if site2 == "":
                site2_dist = "ND"
            else:
                left_site2 = seq_find(site2, start=start0, end=end1)
                if left_site2 == -1:
                    site2_dist = "-1,-1"
                else:
                    right_site2 = seq.rfind(site2, start=start0, end=end1)
                    site2_dist = "%d,%d" % (left_site2 - start0,
                                            end1 - right_site2 - site2_len)
            if uncalled_n <= max_n_freq and frag_len > min_frag_len:
                out_fasta.write(fa_fmt.format(name, seq_id, start0, end1,
                    seq[start0:(start0 + flank_len)],
                    seq[(end1 - flank_len):end1]))
                out_info.write(info_fmt.format(name, seq_id, start0, end1,
                    "included", site2_dist))
            elif uncalled_n > max_n_freq:
                out_info.write(info_fmt.format(name, seq_id, start0, end1,
                    "excluded[N=%.2f]" % uncalled_n, site2_dist))
                logging.debug("Fragment excluded due to Ns:      %s:%d-%d",
                        seq_id, start0, end1)
            else:
                out_info.write(info_fmt.format(name, seq_id, start0, end1,
                    "excluded[L=%d]" % frag_len, site2_dist))
                logging.debug("Fragment excluded short length:      %s:%d-%d (%d)",
                        seq_id, start0, end1, frag_len)

            left_site = right_site
    log_frag_len_histogram(frag_lengths, max(frag_lengths) // 30 )


def log_frag_len_histogram(length_list, binsize):
    counts = numpy.bincount([x // binsize for x in length_list])
    bin_edges = ["[%5d - %5d[" % (i * binsize, (i + 1) * binsize) for i in range(len(counts))]
    bin_edges[-1] = ">%5d         " % (len(counts) * binsize)
    scale = 40.0 / max(counts) 
    logging.info("Size distribution of included fragments (mean = %d):",
            numpy.mean(length_list))
    for i in range(len(counts)):
        s = int(counts[i] * scale) * "*"
        logging.info("%s |%-45s|--(%5d)", bin_edges[i], s, counts[i])

