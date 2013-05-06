__author__ = 'will'

from nose.tools import ok_, eq_, raises
from GeneralSeqTools import fasta_reader
import HIVTransTool
from itertools import product


def testmapping():

    ref_path = 'HIVDBFiles/HXB2Sequence.fasta'
    with open(ref_path) as handle:
        ref_seq = list(fasta_reader(handle))[0][1]

    starts = [1000, 3000, 5000, 6000, 7000]
    widths = [60, 100, 500, 700, 1000]

    check_seqs = []

    ntmp = '%i-%i'
    for start, width in product(starts, widths):
        stop = start+width
        check_seqs.append((ntmp % (start, stop), ref_seq[start:stop]))

    count = 0
    for name, start, seq in HIVTransTool.map_seqs_to_ref(check_seqs, ref_path):
        count += 1
        start = int(start)
        stop = start + len(seq)
        nstart, nstop = map(int, name.split('-'))
        eq_(nstart, start-1, 'Did not find the correct start')
        eq_(nstop, stop-1, 'Did not find the correct end')
        cseq = ref_seq[nstart:nstop]
        eq_(seq, cseq, 'Sequence was wrong')
    eq_(count, len(check_seqs), 'Did not find all of the matches')