__author__ = 'will'
from sklearn.base import BaseEstimator
from sklearn.feature_extraction import DictVectorizer
from sklearn.neighbors import KNeighborsClassifier
from collections import Counter
import numpy as np
from GeneralSeqTools import fasta_reader


class SeqTransformer(BaseEstimator):

    def __init__(self):
        pass

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        max_len = max(len(seq) for seq in X)
        Xout = []
        for row in X:
            Xout.append(row.ljust(max_len, '-'))
        return np.array(Xout)

    def fit_transform(self, X, y=None):
        return self.transform(X)

    @staticmethod
    def get_from_fasta_handle(handle, letters_only=True):

        names = []
        seqs = []
        for name, seq in fasta_reader(handle):
            names.append(name)
            if letters_only:
                seqs.append(''.join(l for l in seq if l.isalpha()))
            else:
                seqs.append(seq)
        return names, SeqTransformer().fit_transform(seqs)


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
        self.keep_mask = None

    def fit(self, X, y=None):
        mask = np.array(map(lambda x: not all(not l.isalpha() for l in x), X.flatten()))
        self.keep_mask = mask.reshape(X.shape)
        return self

    def transform(self, X, y=None):

        if self.keep_mask is None:
            self.fit(X)

        if X.shape[0] != self.keep_mask.shape[0]:
            raise NotImplementedError
        if X.shape[1] != self.keep_mask.shape[1]:
            raise NotImplementedError
        return X[self.keep_mask].reshape(-1, 1)

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X, y)

    def reverse_transform(self, X):

        if np.prod(X.shape) != np.sum(self.keep_mask.flatten()):
            raise NotImplementedError

        Xout = [np.nan]*np.prod(self.keep_mask.shape)
        for pos, val in zip(np.where(self.keep_mask.flatten())[0], X.flatten()):
            Xout[pos] = val

        return np.array(Xout).reshape(self.keep_mask.shape)


class KMerTransform(BaseEstimator):

    def __init__(self, min_k=2, max_k=5, dv=DictVectorizer()):

        self.min_k = min_k
        self.max_k = max_k
        self.dv = dv

    def _generate_kmer(self, seq, k):
        counter = Counter(seq[s:s+k] for s in range(len(seq)-k+1))
        rm_keys = [key for key in counter.keys() if any(not l.isalpha() for l in key)]
        for key in rm_keys:
            counter.pop(key, None)
        return counter

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


class AlignHIV(BaseEstimator):

    def __init__(self, winsize=30, min_k=2, max_k=5):
        self.winsize = winsize
        self.min_k = min_k
        self.max_k = max_k

    def fit(self, X, y):
        """X should be a set of
        """

