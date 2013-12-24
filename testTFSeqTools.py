__author__ = 'will'
from nose.tools import ok_, eq_
import os
import TFSeqTools


def test_load_pwms():

    pwm_dict = TFSeqTools.Load_PWMS()
    all_found = True
    missing = []
    all_correct = True
    wrongs = []
    direc = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(direc, 'HIVDBFiles', 'Jaspar_PWMs.txt')
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name = line.strip().split()[-1].lower()
                if name not in pwm_dict:
                    all_found = False
                    missing.append(name)
                elif name + '-R' not in pwm_dict:
                    all_found = False
                    missing.append(name)
                elif 'motif' not in str(type(pwm_dict[name])):
                    all_correct = False
                    wrongs.append(name)

    yield ok_, all_found, 'Missing: ' + ', '.join(missing)
    yield ok_, all_correct, 'Wrong: ' + ', '.join(wrongs)


def test_reverse_motif():

    pwm_dict = TFSeqTools.Load_PWMS()
    for key, mot in pwm_dict.items():
        if not key.endswith('-R'):
            yield check_rev_mot, key, mot


def check_rev_mot(key, mot):

    rmot = TFSeqTools.true_motif_rev_complement(mot)
    rscore = rmot.pssm.calculate(rmot.consensus)
    fscore = mot.pssm.calculate(mot.consensus)
    eq_(rscore, fscore, '%s was not reversed properly!' % key)


def test_align_to_ref():

    test_base = 'ACTGTTTTGCGTA'
    test_ref = 'ACTGTTgTTGCGTA'

    res_base = 'ACTGTT-TTGCGTA'
    res_ref = 'ACTGTTgTTGCGTA'

    out_base, out_ref = TFSeqTools.align_to_ref(test_ref, test_base)
    tests = [(eq_, res_base.lower(), out_base.lower(),
              'Query alignment was wrong!'),
             (eq_, res_ref.lower(), out_ref.lower(),
              'Reference alignment was wrong!'),
             (ok_, (test_ref, test_base) in TFSeqTools.align_to_ref,
              'Memoization didnt save properly'),
             (eq_, TFSeqTools.align_to_ref[(test_ref, test_base)],
              (out_base, out_ref), 'Memoization saved wrong data!')
             ]
    for tup in tests:
        yield tup


def test_slice_to_ref():

    test_base = 'ACTGTTgggTTGCGTA'
    test_ref = 'ACTGTTTTGCGTA'

    start = 4
    stop = 8

    res = TFSeqTools.slice_to_ref(test_ref, start, stop, test_ref)
    yield eq_, test_ref[4:8].lower(), res.lower(), 'Simple Slice is wrong!'

    expect = 'TTgggTT'
    res = TFSeqTools.slice_to_ref(test_ref, start, stop, test_base)
    yield eq_, expect.lower(), res.lower(), 'Aligned Slice is wrong!'


def test_simple_score_pwm():

    pwm_dict = TFSeqTools.Load_PWMS()
    mot = pwm_dict['arnt']
    #A   4 19  0  0  0  0
    #C  16  0 20  0  0  0
    #G   0  1  0 20  0 20
    #T   0  0  0  0 20  0

    tseq = 'AAACACGTGAAAA'

    cor_seq = 'CACGTG'
    cor_pos = tseq.find(cor_seq)

    _, bpos, nseq = TFSeqTools.simple_score_pwm(mot, tseq)
    yield eq_, bpos, cor_pos, 'Wrong position found!'
    yield eq_, nseq, cor_seq, 'Wrong sequence found!'