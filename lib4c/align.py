"""
Usage:
    4c align [--out=NAME] [--qual=QUAL] <read1> <read2> <config> <index>

Options:
    --out=NAME    output file name [default: stdout]
    --qual=QUAL   fastq quality scale (solexa, phred64, phred33)
                  [default: phred33]

Arguments:
    read1         fastq file of first read in pair
    read2         fastq file of second read in pair
    config        configuration file describing flanks for each bait
    index         bowtie index of restriction enzyme flanks

    Output goes to stdout by default so that any filtering with
    samtools can be done in a pipeline format.
    
Description:
    Take a pair of compressed fastq files, determine which pairs are
    valid (i.e. have the 6-hitter and 4-hitter flank), extract the
    target sequence from the 6-hitter arm (not including the
    restriction site itself), and align the 6-hitter arm against a
    bowtie index containing only sequences flanking the 6-hitter
    sites.  Output is in bam format.

    This alignment takes 16 cores, so it needs to be run on a 
    processing node.
    
Config file format:
    # comment line
    sample name:enzyme1 flank:enzyme2 flank

    enzyme 1: enzyme used for first digest (traditionally a 6 hitter,
        but for higher resolution experiments a 4-hitter).  The fragment
        size of this digest determines the resolution. 
    enzyme 2: enzyme used for second digest
    
    Flank format:
        * flank is separated into 2 parts separated by comma:
            primed_flank,non_primed_flank
          where the non_primed_flank is sequence between the
          primer and the restriction site
        * sequence is upper case except for the RE site, which
          is lower case
        * full sequence of the RE site is included in the flank
          sequence.
"""


import docopt
from . import validators as val
import sys
import os
import shlex
import logging
import subprocess
import collections
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator

BUFSIZE = 81920
#TODO:  refactor such that failure anywhere in the pipeline 
#TODO:      cause all the resources to be shut down properly??
#TODO:      maybe use threading?

def check_args(args):
    schema = {"<read1>": (val.is_valid_infile,
                          "file {<read1>} not found".format(**args)),
              "<read2>": (val.is_valid_infile,
                          "file {<read2>} not found".format(**args)),
              "<config>": (val.is_valid_infile,
                           "config file {<config>} not found".format(**args)),
              "<index>": (val.is_valid_bowtie_index,
                           "{<index>} is not a valid bowtie index".format(**args)),
              "--out": (val.is_valid_outfile,
                        "{--out} is not a valid output file".format(**args)),
              "--qual": (val.is_valid_qscale,
                         "{--qual} is not a valid quality scale".format(**args))}
    ok, errors = val.validate(args, schema)
    if not ok:
        for e in errors:
            logging.error(e)
        return None
    else:
        args["--qual"] = args["--qual"].lower()
        return args

def main(cmdline):
    args = docopt.docopt(__doc__, argv=cmdline)
    args = check_args(args)
    if args is None:
        sys.exit(1)
    else:
        align(args["<read1>"],
              args["<read2>"],
              args["<config>"],
              args["<index>"],
              args["--qual"],
              args["--out"])
    
def align(read1, read2, config, index, qual, out):
    """driver function for the align action"""
    logging.info("***** Starting alignment *****")
    logging.info("Read 1: %s", read1)
    logging.info("Read 2: %s", read2)
    logging.info("Genome: %s", index)
    logging.info("Config: %s", config)
    logging.info("Output: %s", out)
    logging.info("Qual  : --%s-quals", qual)
    try:
        bowtie_version=subprocess.check_output(["bowtie", "--version"])
        logging.info("bowtie:   %s", bowtie_version)
    except OSError:
        logging.error("bowtie executable not found on path")
        sys.exit(1)


    with open(config) as config_fh:
        flanks, flank_prefix_len = parse_config_file(config_fh, min_prefix_len = 6)

    pairs       = make_pairs(read1, read2)
    valid_pairs = process_pairs(pairs, flanks, flank_prefix_len)

    if out == "stdout":
        out_fh = sys.stdout
    else:
        out_fh = open(out, "wb")
    run_bowtie(valid_pairs, index, qual, out_fh)
    out_fh.close()

def run_bowtie(valid_pairs, genome, qual_scale, out_fh):
    """start a bowtie process and feed it's stdin with sequences for
    alignment"""
    bowtie_cmd   = "bowtie --best --all --strata --chunkmbs 256 -m1 --threads=14 --%s-quals --sam %s -" % (
            qual_scale, genome)
    logging.info(bowtie_cmd)
    samtools_cmd = "samtools view -Sb -F4 -"
    logging.info(samtools_cmd)

    bowtie = subprocess.Popen(
            shlex.split(bowtie_cmd),
            bufsize = BUFSIZE,
            stdin   = subprocess.PIPE,
            stdout  = subprocess.PIPE,
            stderr  = subprocess.PIPE)
     
    samtools = subprocess.Popen(
            shlex.split(samtools_cmd),
            bufsize = BUFSIZE,
            stdin   = bowtie.stdout,
            stdout  = out_fh,
            stderr  = subprocess.PIPE,
            close_fds=True)  
    # close_fds is necessary to make sure samtools does not inherit the
    # file descriptors of bowtie which leads to bowtie not terminating
    # see
    # http://stackoverflow.com/questions/1595492/blocks-send-input-to-python-subprocess-pipeline
    
    fq = "@{0}\n{1}\n+\n{2}\n"
    logging.info("START feeding bowtie")
    for pair in valid_pairs:
        try:
            bowtie.stdin.write(fq.format(*pair))
        except IOError:
            logging.error("Bowtie seems to have exited (return code %d)",
                    bowtie.poll())
            logging.error(bowtie.stderr.read())
            samtools.terminate()
            sys.exit(1)
    bowtie.stdin.close()
    logging.info("DONE  feeding bowtie")
    bowtie.wait()
    logging.info("Bowtie return code: %d", bowtie.returncode)
    logging.info(bowtie.stderr.read())
    samtools.wait()
    logging.info("Samtools return code: %d", samtools.returncode)

#TODO: incorporate this class into workflow - it's easier that the tuple i used
#TODO:     until now
class flank(object):
    __slots__ = ("flank6", "flank4", "trim6", "sample")
    def __init__(self, flank6, flank4, trim6, sample):
        self.flank6 = flank6
        self.flank4 = flank4
        self.trim6  = trim6
        self.sample = sample

def _parse_config_file(fh):
    """take a config file and create a list of tuples that contain
    information about the bait flanks for each sample"""
    id_tuples = []
    unique_flanks = set([])
    for line in fh:
        if line.startswith("#"):
            continue
        sample, enz1_flank, enz2_flank = line.strip().split(":")

        enz1_flank_full     = "".join(enz1_flank.split(","))
        enz1_flank_trim_to  = len(enz1_flank_full)
        enz1_flank_re_site  = re.search(r"[gatc]+$", enz1_flank_full)
        if enz1_flank_re_site is not None:
            logging.info("  sample {}: enzyme 1 RE site is {}".format(
                    sample, enz1_flank_re_site.group().upper()))
        else:
            logging.warn("  could not find RE site for enzyme 1 in config file")
        enz1_flank_full     = enz1_flank_full.upper()
        if enz1_flank_full in unique_flanks:
            logging.error("Same enzyme 1 cutter flank occured twice in flank config \
                    file: %s/%s", sample, enz1_flank_full)
            sys.exit(1)
        else:
            unique_flanks.add(enz1_flank_full)

        enz2_flank_full = "".join(enz2_flank.split(",")).upper()
        id_tuples.append((enz1_flank_full,
            enz2_flank_full,
            enz1_flank_trim_to,
            sample))
    return id_tuples

def find_unique_prefix_len(str_lst, min_prefix_len, excl_str_lst):
    """find shortest prefix larger than min_prefix_len that uniquely
    identifies all strings in str_lst; and does not appear in
    excl_str_lst. only works for small-ish data sets"""
    max_word_len = max(len(x) for x in str_lst)
    min_word_len = min(len(x) for x in str_lst)
    if min_word_len < min_prefix_len:
        logging.error("min prefix len is more than shortest flank")
        sys.exit(1)
    for prefix_len in range(min_prefix_len, max_word_len):
        bucket = set([])
        for word in str_lst:
            if len(word) < prefix_len:
                logging.error("flank %s is not unique", word)
                sys.exit(1)
            bucket.add(word[0:prefix_len])
        if len(bucket) == len(str_lst):
            bucket_excl = set([x[0:prefix_len] for x in excl_str_lst])
            if bucket & bucket_excl == set([]):
                return prefix_len
            else:
                continue
    return None

def parse_config_file(fh, min_prefix_len):
    """creates efficient structure for identifying different flanks"""
    flank_tuple = _parse_config_file(fh)
    # length of shortest (but at least min_prefix_len nt long) unique prefix - first enzyme
    # note that the 4-cutter side can be redundant!
    all_enz1_flanks = [x[0] for x in flank_tuple]
    all_enz2_flanks = [x[1] for x in flank_tuple]
    unique_prefix_len = find_unique_prefix_len(all_enz1_flanks, min_prefix_len,
                                               all_enz2_flanks)
    if unique_prefix_len is None:
        logging.error("could not find unique prefix for enzyme1 flanks")
        sys.exit(1)
    else:
        logging.info("enzyme1 unique prefix length: %d", unique_prefix_len)

    # create data structure used for identyfying flanks
    flank_id = {}
    for flank_enz1, flank_enz2, trim_enz1, sample in flank_tuple:
        flank_id[flank_enz1[0:unique_prefix_len]] = (
                flank_enz1, flank_enz2, trim_enz1, sample )
        if len(flank_enz2) > 50:
            logging.debug("%s => (%s, \n     %s..., \n     %d, %s)",
                    flank_enz1[0:unique_prefix_len], flank_enz1, flank_enz2[0:50], trim_enz1, sample)
        else:
            logging.debug("%s => (%s, \n     %s, \n     %d, %s)",
                    flank_enz1[0:unique_prefix_len], flank_enz1, flank_enz2, trim_enz1, sample)
    return flank_id, unique_prefix_len



def parse_gzip_fastq(file_name):
    """uses biopythons FastqGeneralIterator to parse gziped fastq.
    uncompressing is done in a separate process"""
    zcat = subprocess.Popen(["zcat", file_name],
            stdout    = subprocess.PIPE,
            stderr    = subprocess.PIPE,
            bufsize   = BUFSIZE,
            close_fds = True)

    logging.info("START parsing %s", file_name)
    for seq_tuple in FastqGeneralIterator(zcat.stdout):
        yield seq_tuple
    retcode = zcat.wait()
    if retcode != 0:
        logging.error("zcat encountered error (retcode %d)", retcode)
        logging.error(zcat.stderr.read())
        sys.exit(1)
    logging.info("DONE  parsing %s (zcat returned %d)", file_name, retcode)

def make_pairs(read1_file_name, read2_file_name):
    """merge arms of paired end sequences"""
    read1_it = parse_gzip_fastq(read1_file_name)
    read2_it = parse_gzip_fastq(read2_file_name)
    done     = 0
    while True:
        try:
            r1 = read1_it.next()
        except StopIteration:
            done |= 1
        try:
            r2 = read2_it.next()
        except StopIteration:
            done |= 2
        if done == 0:    
            id1 = r1[0].split()[0]
            id2 = r2[0].split()[0]
            if id1[-2] == "/":
                id1 = id1.split("/")[0]
                id2 = id2.split("/")[0]
            if id1 != id2:
                logging.error("Ids of pair do not match: %s vs %s", id1, id2)
                sys.exit(1)
            yield (id1, r1[1], r1[2], r2[1], r2[2])
        elif done == 1:
            logging.error("read1 file contained fewer reads than read2 file")
            sys.exit(1)
        elif done == 2:
            logging.error("read2 file contained fewer reads than read1 file")
            sys.exit(1)
        else:
            break

def has_flank(s, flank, start):
    """determines if string s starts with string flank, allowing up to 3
    mismatches"""
    nMismatch = sum(1 for a, b in zip(s[start:len(flank)], flank[start:]) if a != b)
    if nMismatch <= 3:
        #print "[D] n = %d :%s\n[D]        %s" % (nMismatch, s, flank)
        return True
    else:
        return False

def process_pairs(pair_it, flanks, prefix_len):
    """take a pair iterator, determine which pairs are valid, remove the flanks,
    and yield only the sequence/quality of the enzyme1 target sequence"""
    FLANK_ENZ1, FLANK_ENZ2, TRIM_ENZ1, SAMPLE = range(4)
    total_n, valid_n = 0, 0
    sample_counts = collections.defaultdict(int)
    for _, seq1, qual1, seq2, qual2 in pair_it:
        total_n += 1
        # find which sample the pair matches
        s = flanks.get(seq1[0:prefix_len], None)
        if s is not None:
            # check the two flanks
            if has_flank(seq1, s[FLANK_ENZ1], prefix_len) and \
               has_flank(seq2, s[FLANK_ENZ2], prefix_len):
                valid_n += 1
                new_rid = "%s:%d" % (s[SAMPLE], total_n)
                sample_counts[s[SAMPLE]] += 1
                yield (new_rid, seq1[s[TRIM_ENZ1]:], qual1[s[TRIM_ENZ1]:])
            continue
        s = flanks.get(seq2[0:prefix_len], None)
        if s is not None:
            # check the two flanks
            if has_flank(seq2, s[FLANK_ENZ1], prefix_len) and \
               has_flank(seq1, s[FLANK_ENZ2], prefix_len):
                valid_n += 1
                new_rid = "%s:%d" % (s[SAMPLE], total_n)
                sample_counts[s[SAMPLE]] += 1
                yield (new_rid, seq2[s[TRIM_ENZ1]:], qual2[s[TRIM_ENZ1]:])
    logging.info("Total pairs screened: %8d", total_n)
    logging.info("Valid pairs found:    %8d", valid_n)
    for k, v in sample_counts.items():
        logging.info("  Sample %20s: %8d", k, v)
