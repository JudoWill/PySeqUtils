__author__ = 'will'
from nose.tools import eq_, ok_
import numpy as np
import HIVAlignTools
from sklearn.pipeline import Pipeline
from collections import Counter
from StringIO import StringIO
from GeneralSeqTools import fasta_writer
from itertools import product


def testSeqTransformer():

    inseqs = ['ATGTCG',
              'ATGG',
              'ATGTAHYTD']
    outdata = np.array(['ATGTCG---',
                        'ATGG-----',
                        'ATGTAHYTD'])

    seqformer = HIVAlignTools.SeqTransformer().fit(None)
    out = seqformer.transform(inseqs)
    ok_(np.all(out == outdata))


def testSeqTransformer_from_fasta():

    handle = StringIO()
    inseqs = [('Seq1', 'ATGTCG'),
              ('Seq2', 'ATGG'),
              ('Seq3', 'ATGTAHYTD')]
    fasta_writer(handle, inseqs)
    handle.seek(0)

    outdata = np.array(['ATGTCG---',
                        'ATGG-----',
                        'ATGTAHYTD'])

    names, out = HIVAlignTools.SeqTransformer.get_from_fasta_handle(handle)
    ok_(np.all(out == outdata))
    ok_(all(t == g for t, g in zip(names, ['Seq1', 'Seq2', 'Seq3'])))


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
    eq_(outdata.shape[0], out.shape[0])
    eq_(outdata.shape[1], out.shape[1])
    ok_(np.all(out == outdata))


def testUnrollTransform_withgaps():

    indata = np.array([['ABC', 'BCD', 'CDE'],
                       ['FGH', 'GHI', 'HIJ'],
                       ['KLM', 'LMN', 'MNO'],
                       ['PQR', 'QRS', 'RST'],
                       ['-QR', 'QRS', '---']])
    outdata = np.array(['ABC', 'BCD', 'CDE', 'FGH',
                        'GHI', 'HIJ', 'KLM', 'LMN',
                        'MNO', 'PQR', 'QRS', 'RST',
                        '-QR', 'QRS']).reshape(-1, 1)

    unroller = HIVAlignTools.UnrollTransform()
    out = unroller.fit_transform(indata)
    eq_(outdata.shape[0], out.shape[0])
    eq_(outdata.shape[1], out.shape[1])
    ok_(np.all(out == outdata))


def testUnrollTransform_reverse_tranform():

    indata = np.array([['ABC', 'BCD', 'CDE'],
                       ['FGH', 'GHI', 'HIJ'],
                       ['KLM', 'LMN', 'MNO'],
                       ['-QR', 'QRS', '---'],
                       ['PQR', 'QRS', 'RST']])
    outdata = indata.copy()
    outdata[3, 2] = np.nan

    unroller = HIVAlignTools.UnrollTransform()
    unrolled = unroller.fit_transform(indata)

    rick_rolled = unroller.reverse_transform(unrolled)
    eq_(outdata.shape[0], rick_rolled.shape[0])
    eq_(outdata.shape[1], rick_rolled.shape[1])
    ok_(np.all(outdata == outdata))


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

    eq_(outdata.shape[0], out.shape[0])
    eq_(outdata.shape[1], out.shape[1])
    ok_(np.all(out == outdata))


def testCountKmer():

    indata = 'ABCDEFGABC'
    tups = ['ABC', 'BCD', 'CDE', 'DEF', 'EFG', 'FGA', 'GAB', 'ABC']
    outdata = Counter(tups)

    kmerformer = HIVAlignTools.KMerTransform()
    out = kmerformer._generate_kmer(indata, 3)
    keys = set(outdata.keys()) | set(out.keys())
    for key in keys:
        eq_(outdata[key], out[key])


def testCountKmer_with_gaps():

    indata = 'ABCDEF--GABC'
    tups = ['ABC', 'BCD', 'CDE', 'DEF', 'GAB', 'ABC']
    outdata = Counter(tups)

    kmerformer = HIVAlignTools.KMerTransform()
    out = kmerformer._generate_kmer(indata, 3)
    keys = set(outdata.keys()) | set(out.keys())
    print out
    print outdata
    for key in keys:
        eq_(outdata[key], out[key])


def testKmerTransform():

    indata = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))[:20].reshape(4, 5)
    pipe = Pipeline(steps=[('window', HIVAlignTools.WindowTransformer(winsize=4)),
                           ('unroll', HIVAlignTools.UnrollTransform())])
    out = pipe.transform(indata)
    kmerformer = HIVAlignTools.KMerTransform(min_k=2, max_k=3).fit(out)
    Xout = kmerformer.transform(out)

    eq_(Xout.shape[0], out.shape[0])
    ok_(Xout.shape[1] > out.shape[0])


def test_score_seq():

    seqs = [('atgtag', 'atgtag', 30.0),
            ('atgtag', 'atg', 15),
            ('GTIJ', 'GTIJ', 15),
            ('GTIJAGATS', 'GTIJAGATS', 38)]

    for s1, s2, score in seqs:
        oscore = HIVAlignTools.score_seq(s1, s2)
        yield eq_, score, oscore


def test_score_seqs():

    seqs = [('atgtag', 'atgtag', 30.0),
            ('atgtag', 'atg', 15),
            ('GTIJ', 'GTIJ', 15),
            ('GTIJAGATS', 'GTIJAGATS', 38)]

    former = HIVAlignTools.SeqTransformer()
    known = former.transform(np.array([s for s, _, _ in seqs]))
    guess = former.transform(np.array([s for _, s, _ in seqs]))
    score = sum(v for _, _, v in seqs)

    oscore = HIVAlignTools.score_seqs(known, guess)
    eq_(score, oscore)


def test_get_seq():

    prots = ['env', 'gag', 'nef', 'pol',
             'rev', 'tat', 'vif', 'vpr',
             'vpu']
    typs = ['PRO', 'DNA']
    wanted = list(product(prots, typs))+[('genome', 'DNA'),
                                         ('ltr', 'DNA')]
    for gene, typ in wanted:
        names, seqs = HIVAlignTools.get_seq(gene, typ)
        yield ok_, len(names) > 1, 'Had issue with %s, %s' % (gene, typ)
        yield eq_, len(names), seqs.shape[0], 'Had issue with %s, %s' % (gene, typ)
