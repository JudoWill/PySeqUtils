__author__ = 'will'
from Bio.Seq import Seq
from Bio import Motif
from itertools import groupby
from operator import methodcaller
from StringIO import StringIO
import os


def Load_PWMS(path=None):
    """Loads the JASPAR pwm matrices as a dict."""

    def yield_motifs(path):
        with open(path) as handle:
            for key, lines in groupby(handle, methodcaller('startswith', '>')):
                if key:
                    name = lines.next().strip().split()[-1].lower()
                else:
                    tmp = ''.join(lines)
                    mot = Motif.read(StringIO(tmp), 'jaspar-pfm')
                    yield name, mot
                    yield name+'-R', mot.reverse_complement()

    if path is None:
        direc = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(direc, 'HIVDBFiles', 'Jaspar_PWMs.txt')

    return dict(yield_motifs(path))