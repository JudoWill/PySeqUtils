__author__ = 'will'

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin, ClassifierMixin
from itertools import product
from scipy.sparse import csr_matrix, eye
from sklearn.pipeline import Pipeline
from sklearn.cluster import KMeans, MiniBatchKMeans, AffinityPropagation
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.cross_validation import Bootstrap

from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.base import BaseEstimator, TransformerMixin, ClassifierMixin
from types import IntType, ListType, TupleType
from sklearn.decomposition import PCA, RandomizedPCA, KernelPCA
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest



####Linker functions!


def silhouette_score_linker(predictor, X, y=None):
    clusters = predictor.predict(X)
    if len(set(clusters)) == 1:
        clusters[-1] += 1
    return silhouette_score(X, clusters)


def normalized_mutual_info_score_linker(predictor, X, y):

    clusters = predictor.predict(X)
    return normalized_mutual_info_score(y, clusters)


def normalized_mutual_info_score_scorefunc(X, y):

    scores = []
    pvals = []
    for col in range(X.shape[1]):
        scores.append(normalized_mutual_info_score(X[:,col], y))
        pvals.append(1)

    return np.array(scores), np.array(pvals)


###Tranformers


class BioTransformer(BaseEstimator, TransformerMixin):

    def __init__(self, typ='nuc'):
        self.typ = typ

    def fit(self, *args):
        return self

    def transform(self, X):
        if self.typ == 'nuc':
            letters = 'ACGT-'
        else:
            letters = 'ARNDCEQGHILKMFPSTWYV-'

        nrows, ncols = X.shape
        #out = eye(nrows, ncols*len(letters), format='csr')
        data = []
        rows = []
        cols = []
        for row in range(nrows):
            for num, (col,l) in enumerate(product(range(ncols), letters)):
                if X[row, col].upper()==l:
                    data.append(1)
                    rows.append(row)
                    cols.append(num)

        return csr_matrix((np.array(data), (np.array(rows), np.array(cols))),
                          shape=(nrows, ncols*len(letters)), dtype=float).todense()



### Estimators
class BinBasedCluster(BaseEstimator):

    def __init__(self, bins=[0, 0.5, 1]+range(5, 36)):
        self.bins=bins

    def fit(self, X, y):

        biny = self.bin_data(y)

        self.pred = NearestCentroid().fit(X, biny)
        return self

    def predict(self, X):
        return self.pred.predict(X)

    def score(self, X, y, is_raw=True):
        clusters = self.pred.predict(X)
        if is_raw:
            return adjusted_rand_score(self.bin_data(y), clusters)
        else:
            return adjusted_rand_score(y, clusters)

    def bin_data(self, y):
        return np.digitize(y, self.bins)

    def make_vern_points(self, X, y):

        sel = SelectKBest(score_func=normalized_mutual_info_score_scorefunc)
        sdata = sel.fit_transform(X, y)
        print X.shape, sdata.shape

        pca = PCA(n_components=2)
        pca_trans = pca.fit_transform(sdata)

        biny = self.bin_data(y)

        pred = NearestCentroid().fit(pca_trans, biny)

        x_min, x_max = pca_trans[:, 0].min() - 1, pca_trans[:, 0].max() + 1
        y_min, y_max = pca_trans[:, 1].min() - 1, pca_trans[:, 1].max() + 1
        xx, yy = np.meshgrid(np.linspace(x_min, x_max, 50),
                             np.linspace(y_min, y_max, 50))
        Z = pred.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)

        return pca_trans, biny, xx, yy, Z