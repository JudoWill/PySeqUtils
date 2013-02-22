__author__ = 'will'

from StringIO import StringIO
from nose.tools import eq_, ok_
from pandas import DataFrame
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
    res = list(GeneralSeqTools.call_muscle(seqs))
    eq_(res, aln)


def test_seq_map_to_ref():

    ref_align = 'ATCTCT--ATCT'
    seq_align = 'A-CCCT-AATCT'
    cor_align = 'A-CCCTATCT'

    res = GeneralSeqTools.seq_map_to_ref(seq_align,ref_align)
    eq_(res, cor_align)


def test_seq_align_to_ref():

    ref_seq = 'ATCGATTGC'
    test_seq = 'ATCGATGC'
    cor_mapping = 'ATCGA-TGC'

    inp = [('test1', test_seq)] * 10

    res = list(GeneralSeqTools.seq_align_to_ref(inp, ref_seq))
    result = [('test1', cor_mapping)] * 10

    eq_(res, result)


def test_seq_align_to_ref_multi():

    ref_seq = 'ATCGATTGC'
    test_seq = 'ATCGATGC'
    cor_mapping = 'ATCGA-TGC'

    inp = [('test1', test_seq)] * 10

    res = list(GeneralSeqTools.seq_align_to_ref(inp, ref_seq, max_workers=5))
    result = [('test1', cor_mapping)] * 10

    eq_(res, result)


def test_convert_seqs_to_dataframe():

    indict = {
        'seq1': list('ATCGATTGC'),
        'seq2': list('ATCGATTGC'),
    }
    inseqs = [('seq1', 'ATCGATTGC'), ('seq2', 'ATCGATTGC')]
    tdf = DataFrame(indict).T

    res = GeneralSeqTools.convert_seqs_to_dataframe(inseqs)
    ok_((res == tdf).all().all())


def test_convert_seqDF_to_list():

    indict = {
        'seq1': list('ATCGATTGC'),
        'seq2': list('ATCGATTGC'),
        }
    inseqs = [('seq1', 'ATCGATTGC'), ('seq2', 'ATCGATTGC')]
    tdf = DataFrame(indict).T

    res = GeneralSeqTools.convert_seqDF_to_list(tdf)
    eq_(res, inseqs)
