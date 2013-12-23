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


def test_get_region_simple():

    checks = [('B', 'ltr', 'dna', 'tGGAaGGGcTaaT'),
              ('B', 'env', 'dna', 'ATGagaGtGaagG'),
              ('C', 'env', 'pro', 'MRVRGILRNCQQW'),
              ('B', 'rev', 'dna', 'ATGGCAGgAAGAAGc'),
              ('C', 'rev', 'pro', 'MAGRSGDSDEA'),
              ('B', 'nef', 'dna', 'ATgGgtggcaagt'),
              ('C', 'nef', 'pro', 'MGGKWSKSSIVGW')]

    for sub, reg, alpha, cor_seq in checks:
        check_seq = ConSeqs.GetConSeq(reg, subtype=sub, alphabet=alpha)
        expr = check_seq.lower().startswith(cor_seq.lower())
        msg = 'Error checking Subtype:%s Region:%s Alpha:%s' % (sub, reg, alpha)
        yield ok_, expr, msg


def test_get_region_tough():

    return
    checks = [('B', 'v3', 'pro', 'CTRPNNNTRKSIHIGPGRAFYTTGEIIGDIRQAHC'),
              ('B', 'gp41', 'pro', 'AVGIGAMFLGFLGAAGS')]

    for sub, reg, alpha, cor_seq in checks:
        check_seq = ConSeqs.GetConSeq(reg, subtype=sub, alphabet=alpha)
        expr = check_seq.lower().startswith(cor_seq.lower())
        msg = 'Error checking Subtype:%s Region:%s Alpha:%s' % (sub, reg, alpha)
        ok_(expr, msg + ':'+check_seq[:len(cor_seq)])
        #yield ok_, expr, msg + ':'+check_seq[:len(cor_seq)]