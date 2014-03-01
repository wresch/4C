"""
validators for command line arguments
"""

import os
import sys

def is_valid_infile(fn):
    return os.path.exists(fn) and os.path.isfile(fn)

def is_valid_outfile(fn, allow_existing=False):
    if os.path.exists(fn):
        if os.path.isfile(fn) and allow_existing:
            return True
        else:
            return False
    else:
        pd = os.path.dirname(fn)
        if os.path.exists(pd) and os.path.isdir(pd):
            return True
        else:
            return False
                    
def is_dir(dn, allow_existing=False):
    if os.path.exists(dn) and os.path.isdir(dn):
        if allow_existing:
            return True
        else:
            return False
    else:
        return False
