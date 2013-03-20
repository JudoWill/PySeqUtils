__author__ = 'will'

from nose.tools import ok_
from itertools import product
import Rtools
from pandas import DataFrame, Series, Index
from scipy import randn


def test_quantile_normalize():

    test_data = DataFrame(
        {
            'colA': Series(randn(1, 5000).flatten()),
            'colB': Series(100 * randn(1, 5000).flatten()),
            'colC': Series(100 + randn(1, 5000).flatten()),
        },
    )
    test_data.index = Index(map(str, range(5000)))

    normed_data = Rtools.quantile_norm_with_R(test_data)

    yield ok_, (normed_data.index == test_data.index).all()
    yield ok_, (normed_data.columns == test_data.columns).all()

    means = normed_data.mean(axis=0)

    for a, b in product(means.values, repeat=2):
        yield ok_, abs(a - b) < 0.01


