__author__ = 'will'
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC
from itertools import groupby
from operator import methodcaller
from StringIO import StringIO
import numpy as np
import os
import GeneralSeqTools


class memoize(dict):

    def __init__(self, func, **kwargs):
        super(memoize, self).__init__(**kwargs)
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
                    mot = motifs.read(StringIO(tmp), 'pfm')
                    yield name, mot
                    yield name+'-R', mot.reverse_complement()

    if path is None:
        direc = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(direc, 'HIVDBFiles', 'Jaspar_PWMs.txt')

    return dict(yield_motifs(path))


@memoize
def align_to_ref(ref_seq, base_seq):
    """Aligns a sequence to the reference and caches the result for fast
     lookup later. Returns a tuple (base_seq, ref_seq) properly aligned.

     ref_seq -- The reference sequence to use as a guide.
     query_seq -- The query sequence.

     Returns:
     query_aln -- The aligned query sequence.
     ref_aln -- The aligned reference sequence.
    """

    seqs = [('query', base_seq), ('ref', ref_seq)]
    aligned = dict(GeneralSeqTools.call_muscle(seqs))
    return aligned['query'], aligned['ref']


def slice_to_ref(ref, start, stop, base_seq):
    """Slices the query-sequence based on the reference sequence provided.

     ref_seq -- The reference sequence to use as a guide.
     start, stop -- The start/stop positions to extract. In SLICE syntax!
     query_seq -- The query sequence.

     Returns:
     clipped_seq -- The sequence clipped to the desired region.
    """

    base_align, ref_align = align_to_ref(ref, base_seq)

    ref_pos = 0
    out_seq = ''
    for b_let, r_let in zip(base_align, ref_align):
        ref_pos += r_let != '-'
        if ref_pos > start:
            out_seq += b_let
        if ref_pos == stop:
            return out_seq


def simple_score_pwm(PWM, seq, include_revc=True):
    """Scans a sequence with a PWM and returns the best position, sequence,
     and score.

     PWM -- A Bio.Motif object
     seq -- A sequence to scan
     include_revc=True -- Include a scan with the reverse-complement of the
                          motif.

     Returns:
     score -- The best score for the PWM found in the region.
     pos -- The starting position of the best scoring motif.
     matched_seq --
    """

    bseq = Seq(seq, alphabet=IUPAC.unambiguous_dna)
    scores = PWM.scanPWM(bseq)
    bpos = np.argmax(scores)
    bscore = scores[bpos]

    if include_revc:
        rev_scores = PWM.reverse_complement().scanPWM(bseq)
        if np.max(rev_scores) > bscore:
            bscore = np.max(rev_scores)
            bpos = np.argmax(rev_scores)
    nseq = seq[bpos:(bpos+len(PWM))]
    return bscore, bpos, nseq