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
