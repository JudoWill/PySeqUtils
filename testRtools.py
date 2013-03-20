__author__ = 'will'

from nose.tools import ok_, eq_, raises
from itertools import product
import Rtools
from pandas import DataFrame, Series, Index
from scipy import randn
import pandas.rpy.common as com
import rpy2.robjects as robjects


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


def test_converting_to_factors():

    test_data = DataFrame(
        {
            'colA': Series(randn(1, 5000).flatten() > 0),
            'colB': Series(100 * randn(1, 5000).flatten()),
            'colC': Series(100 + randn(1, 5000).flatten()),
            'colD': Series(randn(1, 5000).flatten() > 0),
        },
    )

    test_data['colA'] = test_data['colA'].map(str)
    test_data['colD'] = test_data['colD'].map(str)

    factor_cols = [('colA', 'True'),
                   ('colD', 'True')]

    rpy_test_df = com.convert_to_r_dataframe(test_data)

    rpy_out_df = Rtools.convert_columns_to_factors(rpy_test_df, factor_cols)
    test_cols = [('colA', 'factor'),
                 ('colB', 'numeric'),
                 ('colC', 'numeric'),
                 ('colD', 'factor')]

    for col, typ in test_cols:
        if typ == 'factor':
            yield eq_, rpy_out_df.rx2(col).nlevels, 2
        elif typ == 'numeric':
            yield ok_, (not hasattr(rpy_out_df.rx2(col), 'nlevels'))

