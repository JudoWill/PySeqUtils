__author__ = 'will'

from nose.tools import ok_, eq_, raises
from GeneralSeqTools import fasta_reader
import HIVTransTool
from itertools import product
import csv

ref_path = 'HIVDBFiles/HXB2Sequence.fasta'


def test_known_mappings():

    with open('TestData/LocatorRes.tsv') as handle:
        cor_res = list(csv.DictReader(handle, delimiter='\t'))

    with open('TestData/testSeqs.fasta') as handle:
        test_seqs = list(fasta_reader(handle))

    for row, crow in zip(HIVTransTool.process_seqs(test_seqs, extract_regions=True), cor_res):
        for f in crow.keys():
            if row[f] is None:
                row[f] = ''
            yield eq_, str(row[f]), crow[f], f


def test_extract_region():

    test_seqs = [
        ('AAAABBBBCCCC', (0, 12), (0, 12), 'AAAABBBBCCCC'),
        ('AAAABBBBCCCC', (0, 12), (4, 8), 'BBBB'),
        ('AAAABBBBCCCC', (4, 16), (3, 7), '-AAA'),
        ('AAAABBBBCCCC', (4, 16), (0, 4), '----'),
        ('AAAABBBBCCCC', (4, 16), (14, 18), 'CC--'),
    ]

    for found_seq, (fs, fe), (rs, re), cor in test_seqs:
        out = HIVTransTool.extract_region(found_seq, fs, fe, rs, re)
        eq_(out, cor, '(%i, %i), (%i, %i)' % (fs, fe, rs, re))