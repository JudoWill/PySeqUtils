__author__ = 'will'
from nose.tools import eq_, ok_
import numpy as np
import HIVAlignTools


def testWindowTransformer():

    winsize = 3
    indata = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))[:20].reshape(4, 5)
    outdata = np.array([['ABC', 'BCD', 'CDE'],
                        ['FGH', 'GHI', 'HIJ'],
                        ['KLM', 'LMN', 'MNO'],
                        ['PQR', 'QRS', 'RST']])

    winformer = HIVAlignTools.WindowTransformer(winsize=winsize)
    winformer.fit(indata)
    out = winformer.transform(indata)
    ok_(np.all(out==outdata))