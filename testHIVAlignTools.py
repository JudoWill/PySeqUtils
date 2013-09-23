__author__ = 'will'
from nose.tools import eq_, ok_
import numpy as np
import HIVAlignTools
from sklearn.pipeline import Pipeline


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
    ok_(np.all(out == outdata))


def testUnrollTransform():

    indata = np.array([['ABC', 'BCD', 'CDE'],
                       ['FGH', 'GHI', 'HIJ'],
                       ['KLM', 'LMN', 'MNO'],
                       ['PQR', 'QRS', 'RST']])
    outdata = indata.reshape(-1, 1)

    unroller = HIVAlignTools.UnrollTransform()
    out = unroller.fit_transform(indata)
    ok_(np.all(out == outdata))


def testPipe():

    winsize = 3
    indata = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))[:20].reshape(4, 5)
    window_data = np.array([['ABC', 'BCD', 'CDE'],
                            ['FGH', 'GHI', 'HIJ'],
                            ['KLM', 'LMN', 'MNO'],
                            ['PQR', 'QRS', 'RST']])
    outdata = window_data.reshape(-1, 1)

    pipe = Pipeline(steps=[('window', HIVAlignTools.WindowTransformer(winsize=winsize)),
                           ('unroll', HIVAlignTools.UnrollTransform())])
    out = pipe.transform(indata)
    ok_(np.all(out == outdata))