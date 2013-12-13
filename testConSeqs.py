__author__ = 'will'
from nose.tools import ok_, eq_
import ConSeqs
from itertools import product


def test_get_region_file():

    regions = ['env', 'gag', 'nef',
               'pol', 'rev', 'tat', 'vif',
               'vpr', 'vpu']
    typs = ['pro', 'dna']

    for reg, typ in list(product(regions, typs))+[('ltr', 'dna')]:
        msg = '%s:%s' % (reg, typ)
        expr = 'HIV1_CON_200' in ConSeqs.get_region_file(reg, typ)
        yield ok_, expr, msg