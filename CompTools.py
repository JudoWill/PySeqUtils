from __future__ import division
__author__ = 'will'
from collections import Counter
from itertools import product, izip
import numpy as np


def identity_score(a, b, weight=0, null=-1):
    """Returns a particular value whenever it encounters identical inputs."""

    return weight if a.lower() == b.lower() else null


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