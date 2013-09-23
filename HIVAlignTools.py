__author__ = 'will'
from sklearn.base import BaseEstimator
from sklearn.feature_extraction import DictVectorizer
from sklearn.neighbors import KNeighborsClassifier
from collections import Counter
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


class KMerTransform(BaseEstimator):

    def __init__(self, min_k=2, max_k=5, dv=DictVectorizer()):

        self.min_k = min_k
        self.max_k = max_k
        self.dv = dv

    def _generate_kmer(self, seq, k):
        return Counter(seq[s:s+k] for s in range(len(seq)-k+1))

    def _yield_kmer(self, X):

        for row in range(X.shape[0]):
            counter = Counter()
            seq = X[row, :][0]
            for k in range(self.min_k, self.max_k+1):
                counter += self._generate_kmer(seq, k)
            yield counter

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):

        return self.dv.fit_transform(list(self._yield_kmer(X)))


class FindInHIV(BaseEstimator):
    """ Just a simple wrapper around the KNeighborsClassifier.
    """

    def __init__(self, estimator=KNeighborsClassifier(n_neighbors=1)):
        self.estimator = estimator

    def fit(self, X, y):
        self.estimator.fit(X, y)
        return self

    def predict(self, X):
        return self.estimator.predict(X)

    def transform(self, X):
        return self.predict(X)

