__author__ = 'will'
import os, os.path
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.feature_extraction import DictVectorizer
from sklearn.neighbors import KNeighborsClassifier
from collections import Counter
import numpy as np
import pandas as pd
from GeneralSeqTools import fasta_reader, fasta_writer
from subprocess import check_call, STDOUT
import shlex
from Bio.Blast import NCBIXML
from StringIO import StringIO
from Bio.Blast.Applications import NcbiblastxCommandline, NcbiblastnCommandline
from tempfile import NamedTemporaryFile
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.grid_search import GridSearchCV
import pickle

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
                tseq = ''.join(l for l in X[row, start:] if l.isalpha())
                if len(tseq) < self.winsize:
                    tseq = tseq.ljust(self.winsize, '-')
                Xout.append(tseq[0:self.winsize])
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


def score_seq(known, guess, gapopen=10, gapextend=1):

    cmd = 'needle -asequence %(cb)s -bsequence %(seq)s -aformat score -gapopen %(go)f -gapextend %(ge)s -outfile %(out)s'
    with NamedTemporaryFile() as conb_handle:
        fasta_writer(conb_handle, [('SeqA', known)])
        conb_handle.flush()
        os.fsync(conb_handle.fileno())
        with NamedTemporaryFile() as seq_handle:
            fasta_writer(seq_handle, [('Seq1', guess)])
            seq_handle.flush()
            os.fsync(seq_handle.fileno())
            with NamedTemporaryFile() as out_handle:
                param_dict = {'cb': conb_handle.name,
                              'seq': seq_handle.name,
                              'out': out_handle.name,
                              'go': gapopen,
                              'ge': gapextend
                              }
                cmd_list = shlex.split(cmd % param_dict)
                check_call(cmd_list, stdout=open(os.devnull, 'w'), stderr=STDOUT)
                for line in out_handle:
                    parts = line.split()
                    if len(parts) == 4:
                        return float(parts[-1][1:-2])


def score_seqs(known_seqs, guess_seqs, gapopen=10, gapextend=1):

    score = 0.0
    for ind in range(known_seqs.shape[0]):
        score += score_seq(known_seqs[ind], guess_seqs[ind],
                           gapopen=gapopen, gapextend=gapextend)
    return score


class BlastAligner(BaseEstimator, ClusterMixin):

    def __init__(self, evalue=10, word_size=2, gapopen=11, gapextend=1,
                 max_intron_length=20, tmp_path='/tmp/', result_type='pro',
                 db_path=NamedTemporaryFile(suffix='.fasta').name, num_threads=1):
        self.evalue = evalue
        self.word_size = word_size
        self.gapopen = gapopen
        self.gapextend = gapextend
        self.max_intron_length = max_intron_length
        self.tmp_path = tmp_path
        self.result_type = result_type
        self.db_path = db_path
        self.num_threads = num_threads

    def _write_seqs(self, X, handle):

        seqs = []
        for row in range(X.shape[0]):
            seq = ''.join(X[row])
            seqs.append(('Seq-%03i' % row, ''.join(l for l in seq if l.isalpha())))

        fasta_writer(handle, seqs)
        handle.flush()
        os.fsync(handle.fileno())

    def fit(self, X, y):
        empty_mask = y == 'XX'

        with open(self.db_path, 'w') as handle:
            self._write_seqs(y[~empty_mask], handle)
        cmd = 'makeblastdb -in %s -dbtype ' % self.db_path
        if self.result_type == 'pro':
            cmd += 'prot'
        else:
            cmd += 'nucl'

        check_call(shlex.split(cmd), stdout=open(os.devnull, 'w'), stderr=STDOUT)
        return self

    def predict(self, X):

        if self.result_type == 'pro':
            blast_cmd = NcbiblastxCommandline
        else:
            blast_cmd = NcbiblastnCommandline

        with NamedTemporaryFile(dir=self.tmp_path, delete=True) as fasta_handle:
            with NamedTemporaryFile(dir=self.tmp_path, delete=True) as out_handle:
                self._write_seqs(X, fasta_handle)
                blastx_cline = blast_cmd(query=fasta_handle.name,
                                         db=self.db_path, outfmt=5,
                                         out=out_handle.name,
                                         evalue=self.evalue,
                                         word_size=self.word_size,
                                         gapopen=self.gapopen,
                                         gapextend=self.gapextend,
                                         max_intron_length=self.max_intron_length,
                                         num_threads=self.num_threads)
                stdout, stderr = blastx_cline()

                blast_records = list(NCBIXML.parse(out_handle))
        prots = []
        for rec in blast_records:
            for align in rec.alignments:
                hsp = align.hsps[0]
                prots.append({'ID': rec.query,
                              'Seq': hsp.query})
        if prots:
            blast_out = pd.DataFrame(prots).groupby('ID')['Seq'].first()
            wanted_out = pd.DataFrame({'ID': ['Seq-%03i' % i for i in range(X.shape[0])],
                                       'want_seq': [True]*X.shape[0],
                                       }).groupby('ID')['want_seq'].first()
            out, _ = blast_out.align(wanted_out, join='right')

            return SeqTransformer().transform(out.fillna('XX').values)
        else:
            return np.array(['XX']*X.shape[0])

    def score(self, X, y):

        empty_mask = y == 'XX'
        out_aligns = self.predict(X)

        pos_scores = score_seqs(y[~empty_mask], out_aligns[~empty_mask])
        bad_scores = score_seqs(out_aligns[empty_mask], out_aligns[empty_mask])
        return (pos_scores - bad_scores)/y.shape[0]


def get_seq(prot_name, typ):
    trans_path = '/home/will/PySeqUtils/TransToolStuff/'
    tmp = 'HIV1_ALL_2012_%s_%s.fasta' % (prot_name.lower(), typ.upper())
    with open(trans_path + tmp) as handle:
        return SeqTransformer.get_from_fasta_handle(handle)


def generate_traindata(prot, train_type='pro'):

    neg_controls = {'env': ['gag', 'pol', 'vif', 'vpr', 'ltr'],
                    'gag': ['ltr', 'vif', 'vpr', 'vpu', 'tat', 'rev', 'env'],
                    'ltr': ['gag', 'pol', 'vpr', 'vpu', 'env'],
                    'nef': ['pol', 'gag', 'vpu', 'tat'],
                    'pol': ['env', 'vpr', 'vpu', 'nef', 'rev', 'ltr'],
                    'rev': ['ltr', 'gag', 'pol', 'vif', 'nef'],
                    'tat': ['ltr', 'pol', 'vif', 'nef'],
                    'vif': ['ltr', 'tat', 'vpu', 'rev', 'env', 'nef'],
                    'vpr': ['ltr', 'gag', 'pol', 'rev', 'env', 'nef'],
                    'vpu': ['ltr', 'gag', 'pol', 'vif', 'vpr', 'nef'],
                    'gp41': ['gag', 'pol', 'vif', 'vpr', 'ltr'],
                    'gp120': ['gag', 'pol', 'vif', 'vpr', 'ltr'],
                    'v3': ['gag', 'pol', 'vif', 'vpr', 'ltr'],
                    }

    pos_names, pos_X = get_seq('genome', 'DNA')
    env_names, env_y = get_seq(prot, train_type)
    neg_names = []
    neg_X = None
    for neg_prot in neg_controls[prot]:
        tnames, tx = get_seq(neg_prot, 'DNA')
        neg_names += tnames
        if neg_X is None:
            neg_X = tx.copy()
        else:
            neg_X = np.concatenate((neg_X, tx))

    pos_X_ser = pd.Series(pos_X, index=pos_names)
    env_y_ser = pd.Series(env_y, index=env_names)

    X_ser, y_ser = pos_X_ser.align(env_y_ser, join='inner')
    X = X_ser.values
    y = y_ser.values
    in_env = set(env_names)
    neg_inds = [num for num, name in enumerate(neg_names) if name not in in_env]
    wneg_X = neg_X[neg_inds]
    wneg_y = np.array(['XX']*wneg_X.shape[0])

    return np.concatenate((X, wneg_X)), np.concatenate((y, wneg_y))


def train_aligner(prot, path, train_type='pro', test_size=500, n_jobs=1, verbose=0):

    X, y = generate_traindata(prot, train_type=train_type)
    yclass = y == 'XX'
    cv = StratifiedShuffleSplit(yclass,
                                n_iter=3,
                                test_size=test_size,
                                train_size=100)
    param_dict = {'evalue': np.logspace(-100, 1, 20)}
    gd = GridSearchCV(BlastAligner(num_threads=20,
                                   result_type=train_type),
                      param_dict,
                      cv=cv,
                      n_jobs=n_jobs,
                      verbose=verbose,
                      refit=False)
    gd.fit(X, y)

    aligner = BlastAligner(num_threads=100,
                           result_type=train_type,
                           db_path=path, **gd.best_params_)
    aligner.fit(X, y)

    return aligner


def build_aligners(base_path='/home/will/PySeqUtils/TransToolStuff/dbs/',
                   verbose=1, n_jobs=1):

    prots = ['gag', 'vif',  'vpu',
             'tat', 'rev', 'env', 'gp120',
             'gp41', 'v3', 'ltr'] #'vpr',
    for prot in prots:
        print prot
        path = base_path+prot
        if os.path.exists(path + '.pkl'):
            continue
        aligner = train_aligner(prot, path+'.fasta',
                                verbose=verbose, n_jobs=n_jobs,
                                train_type='pro' if prot != 'ltr' else 'dna')
        with open(path+'.pkl', 'w') as handle:
            pickle.dump(aligner, handle)

if __name__ == '__main__':
    pass