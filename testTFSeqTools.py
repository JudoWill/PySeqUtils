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
                elif 'Bio.Motif._Motif.Motif' not in str(type(pwm_dict[name])):
                    all_correct = False
                    wrongs.append(name)

    yield ok_, all_found, 'Missing: ' + ', '.join(missing)
    yield ok_, all_correct, 'Wrong: ' + ', '.join(wrongs)


def test_align_to_ref():

    test_base = 'ACTGTTTTGCGTA'
    test_ref = 'ACTGTTgTTGCGTA'

    res_base = 'ACTGTT-TTGCGTA'
    res_ref = 'ACTGTTgTTGCGTA'

    out_base, out_ref = TFSeqTools.align_to_ref(test_base, test_ref)
    yield eq_, res_base.lower(), out_base.lower()
    yield eq_, res_ref.lower(), out_ref.lower()
    yield ok_, (test_base, test_ref) in TFSeqTools.align_to_ref


def test_slice_to_ref():

    test_base = 'ACTGTTgggTTGCGTA'
    test_ref = 'ACTGTTTTGCGTA'

    start = 4
    stop = 8

    res = TFSeqTools.slice_to_ref(test_ref, test_ref, start, stop)
    yield eq_, test_ref[4:8].lower(), res.lower(), 'Simple Slice is wrong!'

    expect = 'TTgggTT'
    res = TFSeqTools.slice_to_ref(test_base, test_ref, start, stop)
    yield eq_, expect.lower(), res.lower(), 'Aligned Slice is wrong!'