from __future__ import division
__author__ = 'will'
from collections import Counter
from itertools import product, izip, groupby
from StringIO import StringIO
from operator import methodcaller
import re
import os
import numpy as np


def identity_score(a, b, weight=0, null=-1):
    """Returns a particular value whenever it encounters identical inputs."""

    return weight if a.upper() == b.upper() else null


def replacement_mat_score(dist_mat, a, b, missing=None):
    """Calculates the replacement score based on the provided dict.
    Useful to use with functools.partial!
    """

    return dist_mat.get((a.upper(), b.upper()), missing)


def null_score(gAcol, gBcol, score_func=identity_score):
    """Calculates the expected score under the null-hypothesis that the two
     inputs come from the same distribution.

     It does this by first grouping both sets into one group and then doing a
     pairwise comparison using the identity function.
    """

    counts = Counter(gAcol)+Counter(gBcol)

    scores = []
    weights = []
    for (kA, vA), (kB, vB) in product(counts.items(), repeat=2):
        scores.append(score_func(kA, kB))
        weights.append(vA*vB)

    return np.average(scores, weights=weights)


def group_score(gAcol, gBcol, score_func=identity_score):
    """Calculates the score under the hypothesis that the two inputs come from
     different distributions.

     It does this by doing a pairwise comparison using the identity function.
    """

    pairwise_items = product(gAcol, gBcol)
    pair_counts = Counter(pairwise_items)
    order = sorted(pair_counts.keys())
    nums = np.array([score_func(a, b) for a, b in order])
    weights = np.array([pair_counts[key] for key in order])
    mu = np.average(nums, weights=weights)

    return mu


def group_score_seq(groupA, groupB, score_func=identity_score, has_names=True):
    """Takes an input of multi-aligned sequences and evaluates the group_score
     and null_score at each column in the alignment.
    """

    if has_names:
        Aseqs = [seq for _, seq in groupA]
        Bseqs = [seq for _, seq in groupB]
    else:
        Aseqs = list(groupA)
        Bseqs = list(groupB)

    seq_len = len(Aseqs[0])
    assert all(len(seq) == seq_len for seq in Aseqs), 'All sequences need to be the same length!'
    assert all(len(seq) == seq_len for seq in Bseqs), 'All sequences need to be the same length!'

    Acol_wise = izip(*Aseqs)
    Bcol_wise = izip(*Bseqs)
    mus = []
    nmus = []
    for gAcol, gBcol in izip(Acol_wise, Bcol_wise):
        mus.append(group_score(gAcol, gBcol, score_func=score_func))
        nmus.append(null_score(gAcol, gBcol, score_func=score_func))

    return np.array(mus), np.array(nmus)


def flatten_mat(handle):
    """A simple utility function for flattening the substitution matrix inputs.
    """

    for row_pos, line in enumerate(handle):
        parts = line.split()
        for col_pos, num in enumerate(parts):
            yield row_pos, col_pos, num


def load_sub_mat(instr):
    """Just a testing sub right now
    """

    grab_info = lambda x: x.strip().split(None, 1)[-1]

    ID, desc, ref, author, text, row_line = [None]*6
    handle = StringIO(instr)
    for line in handle:
        if line.startswith('H '):
            ID = grab_info(line)
        elif line.startswith('D '):
            desc = grab_info(line)
        elif line.startswith('R '):
            ref = grab_info(line)
        elif line.startswith('A '):
            author = grab_info(line)
        elif line.startswith('T '):
            text = grab_info(line)
        elif line.startswith('M '):
            row_line = line
            break

    rows, cols = re.findall('= ([\w\-]+)', row_line)
    out_dict = {}
    for row_pos, col_pos, num in flatten_mat(handle):
        out_dict[(rows[row_pos], cols[col_pos])] = float(num)
        out_dict[(cols[col_pos], rows[row_pos])] = float(num)

    return ID, desc, ref, author, text, out_dict


def load_dist_mats(path=None):

    if path is None:
        direc = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(direc, 'HIVDBFiles', 'aa_sub_mat.txt')

    with open(path) as handle:
        keyfunc = methodcaller('startswith', '//')
        for key, lines in groupby(handle, key=keyfunc):
            if not key:
                instr = ''.join(lines)
                yield load_sub_mat(instr)