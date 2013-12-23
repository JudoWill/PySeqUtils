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


def get_region_span(region, alphabet):

    if (alphabet == generic_dna) or (alphabet.lower() == 'dna'):
        region_dict = {
            'gp120': ('env', 0, 7758-6225),
            'gp41': ('env', 7758-6225, 8795-6225),
            'v1': ('env', 393, 471),
            'v2': ('env', 471, 588),
            'v3': ('env', 267*3, 302*3),
            'v4': ('env', 1155, 1254),
            'v5': ('env', 1380, 1407),
            'tat-1': ('tat', 0, 6045-5831),
            'tat-2': ('tat', 6045-5831, 8469-8379),
            'protease': ('pol', 2253-2085, 2550-2085),
            'rt': ('pol', 2550-2085, 3870-2085),
            'rnase': ('pol', 3870-2085, 4230-2085),
            'integrase': ('pol', 4230-2085, 5096-2085),
        }
        return region_dict[region.lower()]
    elif (alphabet == generic_protein) or (alphabet.lower() == 'pro'):
        reg, start, stop = get_region_span(region, 'dna')
        return reg, start/3, stop/3


def GetConSeq(region, alphabet=generic_dna, subtype='B', drop_gaps=True):

    if (alphabet == generic_dna) or (alphabet.lower() == 'dna'):
        path = get_region_file(region, 'dna')
    elif (alphabet == generic_protein) or (alphabet.lower() == 'pro'):
        path = get_region_file(region, 'pro')
    else:
        raise(TypeError, 'alphabet must be: "dna", "pro", generic_dna, generic_protein')
    seq = None
    if path is None:
        new_region, start, stop = get_region_span(region, alphabet)
        conB_seq = GetConSeq(new_region,
                             subtype='B',
                             alphabet=alphabet,
                             drop_gaps=False)
        sub_seq = GetConSeq(new_region,
                            subtype=subtype,
                            alphabet=alphabet,
                            drop_gaps=False)
        nstart, nstop = (None, None)
        print conB_seq
        print sub_seq
        conb_pos = 0
        for aln_pos, l in enumerate(conB_seq):
            if l != '-':
                conb_pos += 1
            if conb_pos == start:
                nstart = aln_pos
            if conb_pos == stop:
                nstop = aln_pos
                break
        seq = sub_seq[nstart:nstop]
    else:

        wanted_key = 'CONSENSUS_'+subtype
        with open(path) as handle:
            for name, seq in GeneralSeqTools.fasta_reader(handle):
                name = name.split('(')[0]
                if name == wanted_key:
                    break

    if drop_gaps:
        return seq.replace('-', '').replace('$', '')
    else:
        return seq.replace('$', '')