__author__ = 'will'

import os
import glob

FILEPATH = os.path.abspath(__file__).rsplit(os.sep, 1)[0]
CONPATH = os.path.join(FILEPATH, 'HIVDBFiles', 'ConSeqFiles')


def get_region_file(region, typ):
    """Returns the path to relevant file."""

    con_files = glob.glob(os.path.join(CONPATH, '*.fasta'))
    for f in con_files:
        if (region.lower() in f) and (typ.upper() in f):
            return f

def GetConSeq(region, subtype='B'):

    pass