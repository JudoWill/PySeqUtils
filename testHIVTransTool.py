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

    for row, crow in zip(HIVTransTool.map_seqs_to_ref(test_seqs), cor_res):
        for f in row.keys():
            if row[f] is None:
                row[f] = ''
            yield eq_, str(row[f]), crow[f], f
