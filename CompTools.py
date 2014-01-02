from __future__ import division
__author__ = 'will'
from collections import Counter
from itertools import product
import numpy as np


def identity_score(a, b, weight=1, null=0):
    """Returns a particular value whenever it encounters identical inputs."""

    return weight if a.lower() == b.lower() else null


def expect_score(gAcol, gBcol, score_func=identity_score):
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