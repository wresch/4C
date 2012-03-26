import logging
import os
import sys
import numpy
import subprocess
from Bio import SeqIO

def make_index(args):
    logging.info("***** Creating new genome index *****")
    logging.info("Genome:     %s", args.genome.name)
    logging.info("RE site:    %s", args.site.upper())
    logging.info("RE name:    %s", args.name)
    logging.info("flank:      %d", args.flank)
    
    out_dir = args.name
    if os.path.exists(out_dir):
        logging.error("Target directory %s already exists", out_dir)
        sys.exit(1)
    os.mkdir(out_dir)
    out_fasta_fh = open(os.path.join(out_dir, "%s.fa" % args.name), "w")
    out_info_fh  = open(os.path.join(out_dir, "%s.info" % args.name), "w")

    make_fasta_file(args.genome, args.site, args.name, args.flank, 
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

