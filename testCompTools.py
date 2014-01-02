__author__ = 'will'
from nose.tools import ok_, eq_
import numpy as np
import CompTools


def test_identity_score():

    cis = CompTools.identity_score
    yield eq_, cis('A', 'A'), 1, 'Missed basic IDENT!'
    yield eq_, cis('A', 'B'), 0, 'Missed basic NON-IDENT!'
    yield eq_, cis('A', 'a'), 1, 'Missed basic IDENT! with cases'
    yield eq_, cis('A', 'A', weight=10), 10, 'Missed basic IDENT with weight'
    yield eq_, cis('A', 'B', null=-5), -5, 'Missed basic IDENT with weight'


def test_expect_score():

    colA = 'AAAAT'
    colB = 'TCCCA'

    score = CompTools.null_score(colA, colB)
    cor_score = 0.38  # by hand
    eq_(score, cor_score, 'Wrong result compared to the hand-calculation')


def test_score_groups():

    colA = 'AAAAT'
    colB = 'TCCCA'

    score = CompTools.group_score(colA, colB)
    cor_score = 0.2  # by hand
    eq_(score, cor_score, 'Wrong result compared to the hand-calculation')


def test_group_score_seq():

    groupA = [('name1', 'AA'),
              ('name1', 'AA'),
              ('name1', 'AA'),
              ('name1', 'AA'),
              ('name1', 'TT')]
    groupB = [('name1', 'TT'),
              ('name1', 'CC'),
              ('name1', 'CC'),
              ('name1', 'CC'),
              ('name1', 'AA')]

    mu, nmu = CompTools.group_score_seq(groupA, groupB)

    ok_(np.all(mu == 0.2), 'Group scores are incorrect!')
    ok_(np.all(nmu == 0.38), 'Null scores are incorrect!')

