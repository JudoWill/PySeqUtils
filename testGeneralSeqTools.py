__author__ = 'will'

from StringIO import StringIO

from nose.tools import eq_

import GeneralSeqTools


def test_fasta_reader():

    input_items = ['>test1', 'ATCTGCTAGTCGA', 'ATCGAGTAGT', '>test2', 'ATCGATGC']

    input_seq = '\n'.join(input_items)

    res = list(GeneralSeqTools.fasta_reader(StringIO(input_seq)))
    eq_(len(res), 2)
    eq_(res[0][0], 'test1')
    eq_(res[0][1], 'ATCTGCTAGTCGAATCGAGTAGT')
    eq_(res[1][0], 'test2')
    eq_(res[1][1], 'ATCGATGC')


def test_fasta_writer():

    items = ['>test1', 'ATCTGCTAGTCGAATCGAGTAGT', '>test2', 'ATCGATGC']
    test_seq = '\n'.join(items) + '\n'

    handle = StringIO()
    GeneralSeqTools.fasta_writer(handle, [('test1', 'ATCTGCTAGTCGAATCGAGTAGT'),
                                            ('test2', 'ATCGATGC')])
    handle.seek(0)
    data = handle.read()
    eq_(test_seq, data)


def test_muscle_basic_call():

    seqs = [('test1', 'ATCGATTGC'), ('test2', 'ATCGATGC')]
    aln = [('test1', 'ATCGATTGC'), ('test2', 'ATCGA-TGC')]
    #muscle_res = '\n'.join(['>test1', 'ATCGATTGC', '>test2', 'ATCGA-TGC'])

    res = list(GeneralSeqTools.call_muscle(seqs))
    eq_(res, aln)



