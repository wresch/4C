"""
validators for command line arguments
"""

import os
import sys
import re

def validate(d, schema):
    errors = []
    ok = True
    for key in schema:
        validator, error = schema[key]
        if not (key in d and validator(d[key])):
            errors.append("{}: {}".format(key, error))
            ok = False
    return ok, errors

def is_one_word(s):
    return re.match("^\w+$", s) is not None

def is_seq(s):
    return re.match("^[AaCcGgTt]+$", s) is not None

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def is_valid_qscale(s):
    if s.lower() in ("phred33", "phred64", "solexa"):
        return True
    else:
        return False

def is_valid_infile(fn):
    return os.path.exists(fn) and os.path.isfile(fn)

def is_valid_infile_list(lst):
    for fn in lst:
        if not is_valid_infile(fn):
            return False
    return True

def is_valid_fa_file(fn):
    return is_valid_infile(fn) and os.path.splitext(fn)[1] == ".fa"

def is_valid_bam_file(fn):
    return is_valid_infile(fn) and os.path.splitext(fn)[1] == ".bam"

def is_valid_bam_file_list(l):
    for bf in l:
        if not is_valid_bam_file(bf):
            return False
    return True

def is_valid_bowtie_index(fn):
    fn1 = "{}.1.ebwt".format(fn)
    return is_valid_infile(fn1)


def is_valid_outfile(fn, allow_existing=False):
    if fn == "stdout":
        return True
    if os.path.exists(fn):
        if os.path.isfile(fn) and allow_existing:
            return True
        else:
            return False
    else:
        pd = os.path.dirname(fn)
        if pd == "":
            return True
        elif os.path.exists(pd) and os.path.isdir(pd):
            return True
        else:
            return False
                    
def is_dir(dn, allow_existing=False):
    """returns true if (1) it's a directory and allow_exisiting
    is true or (2) it does not exists (i.e. can be created)"""
    if os.path.exists(dn) and os.path.isdir(dn):
        if allow_existing:
            return True
        else:
            return False
    else:
        return True
