__author__ = 'will'
from nose.tools import ok_, eq_
import CompTools


def test_identity_score():

    cis = CompTools.identity_score
    yield eq_, cis('A', 'A'), 1, 'Missed basic IDENT!'
    yield eq_, cis('A', 'B'), 0, 'Missed basic NON-IDENT!'
    yield eq_, cis('A', 'a'), 1, 'Missed basic IDENT! with cases'
    yield eq_, cis('A', 'A', weight=10), 10, 'Missed basic IDENT with weight'
    yield eq_, cis('A', 'B', null=-5), -5, 'Missed basic IDENT with weight'