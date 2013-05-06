__author__ = 'will'
from nose.tools import ok_, eq_, raises

import QuasiTools


def test_sam_read():

    queries = [
        'H348453.1',
        'H348451.1',
        'H348449.1',
        'H348447.1',
        'H348445.1',
        'H348443.1',
        'H348441.1',
        'H348439.1']

    start_seq = 'CCAGTA'

    with open('TestData/samfile.sam') as handle:
        for row, sname in zip(QuasiTools.SAMreader(handle), queries):
            eq_(row['QNAME'], sname)
            ok_(row['SEQ'].startswith(start_seq))
