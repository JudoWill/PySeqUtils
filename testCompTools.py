__author__ = 'will'
from nose.tools import ok_, eq_
import numpy as np
import CompTools
import os
from itertools import product
from functools import partial


def test_identity_score():

    cis = CompTools.identity_score
    yield eq_, cis('A', 'A'), 0, 'Missed basic IDENT!'
    yield eq_, cis('C', 'A'), -1, 'Missed basic NON-IDENT!'
    yield eq_, cis('A', 'a'), 0, 'Missed basic IDENT! with cases'
    yield eq_, cis('A', 'A', weight=10), 10, 'Missed basic IDENT with weight'
    yield eq_, cis('A', 'B', null=-5), -5, 'Missed basic IDENT with weight'


def test_replacement_mat_score():

    dmat = {('A', 'A'): 0,
            ('C', 'A'): -1}
    cis = partial(CompTools.replacement_mat_score, dmat)
    yield eq_, cis('A', 'A'), 0, 'Missed basic IDENT!'
    yield eq_, cis('C', 'A'), -1, 'Missed basic NON-IDENT!'
    yield eq_, cis('A', 'a'), 0, 'Missed basic IDENT! with cases'
    yield eq_, cis('G', 'A'), None, 'Missed basic IDENT with weight'
    yield eq_, cis('G', 'A', missing=-5), -5, 'Missed basic IDENT with weight'


def test_expect_score():

    colA = 'AAAAT'
    colB = 'TCCCA'

    score, _ = CompTools.null_score(colA, colB)
    cor_score = -0.62  # by hand
    eq_(score, cor_score, 'Wrong result compared to the hand-calculation')


def test_score_groups():

    colA = 'AAAAT'
    colB = 'TCCCA'

    score, _ = CompTools.group_score(colA, colB)
    cor_score = -0.8  # by hand
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

    df = CompTools.group_score_seq(groupA, groupB, has_names=True)

    yield ok_, np.all(df['InterMu'] == -0.8), 'Group scores are incorrect!'
    yield ok_, np.all(df['NullMu'] == -0.62), 'Null scores are incorrect!'

    groupA = [seq for _, seq in groupA]
    groupB = [seq for _, seq in groupB]

    df = CompTools.group_score_seq(groupA, groupB, has_names=False)
    yield ok_, np.all(df['InterMu'] == -0.8), 'Group scores are incorrect!'
    yield ok_, np.all(df['NullMu'] == -0.62), 'Null scores are incorrect!'


def test_load_sub_mat():

    tmp = """
H WILL-TEST
D A test matrix
R PMID: 2051488
A Dampier, W
T A matrix to test my ability to load in the data
J Will dampier's code
M rows = ACGT, cols = ACGT
      0
     -3.      0
     -7.     -1.      0
     -4.      3.      2.      0
      """

    ID, desc, pmid, author, text, out_dict = CompTools.load_sub_mat(tmp)

    cor_ID = "WILL-TEST"
    cor_desc = "A test matrix"
    cor_pmid = "PMID: 2051488"
    cor_author = "Dampier, W"
    cor_text = "A matrix to test my ability to load in the data"
    cor_dict = {('A', 'A'): 0,
                ('A', 'C'): -3,
                ('A', 'G'): -7,
                ('A', 'T'): -4,
                ('C', 'A'): -3,
                ('C', 'C'): 0,
                ('C', 'G'): -1,
                ('C', 'T'): 3,
                ('G', 'A'): -7,
                ('G', 'C'): -1,
                ('G', 'G'): 0,
                ('G', 'T'): 2,
                ('T', 'A'): -4,
                ('T', 'C'): 3,
                ('T', 'G'): 2,
                ('T', 'T'): 0}

    yield eq_, cor_ID, ID
    yield eq_, cor_desc, desc
    yield eq_, cor_pmid, pmid
    yield eq_, cor_author, author
    yield eq_, cor_text, text
    for key, val in cor_dict.items():
        yield eq_, val, out_dict.get(key, None), 'Incorrect value at (%s, %s)' % key


def test_load_dist_mats():

    direc = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(direc, 'HIVDBFiles', 'aa_sub_mat.txt')

    wanted_IDs = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('H '):
                wanted_IDs.append(line.strip().split(None, 1)[-1])

    test_iter = CompTools.load_dist_mats()
    for out_tup, wanted_id in zip(test_iter, wanted_IDs):
        ID, desc, pmid, author, text, out_dict = out_tup
        yield eq_, ID, wanted_id, 'Incorrect ID found'
        #eq_(ID, wanted_id, 'Incorrect ID found')

        missing = []
        for tup in product('ARNDCQEGHILKMFPSTWYV', 'ARNDCQEGHILKMFPSTWYV'):
            if tup not in out_dict:
                missing.append(tup)
        #ok_(len(missing) == 0, '%s was missing %i items!' % (ID, len(missing)))
        yield ok_, len(missing) == 0, '%s was missing %i items!' % (ID, len(missing))
