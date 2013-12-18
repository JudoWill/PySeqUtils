__author__ = 'will'
from Bio.Seq import Seq
from Bio import Motif
from itertools import groupby
from operator import methodcaller
from StringIO import StringIO
import os
import GeneralSeqTools


class memoize(dict):

    def __init__(self, func):
        self.func = func

    def __call__(self, *args):
        return self[args]

    def __missing__(self, key):
        result = self[key] = self.func(*key)
        return result


def Load_PWMS(path=None):
    """Loads the JASPAR pwm matrices as a dict."""

    def yield_motifs(path):
        with open(path) as handle:
            for key, lines in groupby(handle, methodcaller('startswith', '>')):
                if key:
                    name = lines.next().strip().split()[-1].lower()
                else:
                    tmp = ''.join(lines)
                    mot = Motif.read(StringIO(tmp), 'jaspar-pfm')
                    yield name, mot
                    yield name+'-R', mot.reverse_complement()

    if path is None:
        direc = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(direc, 'HIVDBFiles', 'Jaspar_PWMs.txt')

    return dict(yield_motifs(path))


@memoize
def align_to_ref(base_seq, ref_seq):
    """Aligns a sequence to the reference and caches the result for fast
     lookup later. Returns a tuple (base_seq, ref_seq) properly aligned.
    """

    seqs = [('query', base_seq), ('ref', ref_seq)]
    aligned = dict(GeneralSeqTools.call_muscle(seqs))
    return aligned['query'], aligned['ref']


def slice_to_ref(base_seq, ref, start, stop):

    base_align, ref_align = align_to_ref(base_seq, ref)

    ref_pos = 0
    out_seq = ''
    for b_let, r_let in zip(base_align, ref_align):
        ref_pos += r_let != '-'
        if ref_pos > start:
            out_seq += b_let
        if ref_pos == stop:
            return out_seq

