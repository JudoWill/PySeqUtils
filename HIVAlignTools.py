__author__ = 'will'
from sklearn.base import BaseEstimator
import numpy as np


class WindowTransformer(BaseEstimator):
    """ Converts arrays of items into windows of arrays. Currently only
     works on Char-arrays. For example:

     indata = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))[:20].reshape(4, 5)
     outdata = np.array([['ABC', 'BCD', 'CDE'],
                         ['FGH', 'GHI', 'HIJ'],
                         ['KLM', 'LMN', 'MNO'],
                         ['PQR', 'QRS', 'RST']])
    """

    def __init__(self, winsize=30):
        self.winsize = winsize

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):

        Xout = []
        for row in range(X.shape[0]):
            for start in range((X.shape[1]-self.winsize)+1):
                Xout.append(''.join(X[row, start:start+self.winsize]))
        return np.array(Xout).reshape(X.shape[0], -1)


class UnrollTransform(BaseEstimator):
    """ Simply reshapes the input. This is useful for after the WindowTransformer.
    For example:
    indata = np.array([['ABC', 'BCD', 'CDE'],
                       ['FGH', 'GHI', 'HIJ'],
                       ['KLM', 'LMN', 'MNO'],
                       ['PQR', 'QRS', 'RST']])
    outdata = np.array(['ABC',
                         'BCD',
                         'CDE',
                         'FGH',
                         'GHI',
                         ...])
    """

    def __init__(self, numcols=None):
        self.numcols = numcols

    def fit(self, X, y=None):
        self.numcols = X.shape[1]
        return self

    def transform(self, X, y=None):
        return X.reshape(-1, 1)

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X, y)