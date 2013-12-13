__author__ = 'will'

import os
import glob
from Bio.Alphabet import generic_dna, generic_protein

import GeneralSeqTools

FILEPATH = os.path.abspath(__file__).rsplit(os.sep, 1)[0]
CONPATH = os.path.join(FILEPATH, 'HIVDBFiles', 'ConSeqFiles')


def get_region_file(region, typ):
    """Returns the path to relevant file."""

    con_files = glob.glob(os.path.join(CONPATH, '*.fasta'))
    for f in con_files:
        if (region.lower() in f) and (typ.upper() in f):
            return f


def GetConSeq(region, alphabet=generic_dna, subtype='B'):

    if (alphabet == generic_dna) or (alphabet.lower() == 'dna'):
        path = get_region_file(region, 'dna')
    elif (alphabet == generic_protein) or (alphabet.lower() == 'pro'):
        path = get_region_file(region, 'pro')
    else:
        raise(TypeError, 'alphabet must be: "dna", "pro", generic_dna, generic_protein')

    wanted_key = 'CONSENSUS_'+subtype
    with open(path) as handle:
        for name, seq in GeneralSeqTools.fasta_reader(handle):
            name = name.split('(')[0]
            if name == wanted_key:
                return seq.replace('-', '').replace('$', '')