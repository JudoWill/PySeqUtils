__author__ = 'will'
import os, os.path
import csv
import shlex
from itertools import ifilter
from subprocess import check_output, check_call
from GeneralSeqTools import fasta_writer
from TreeingTools import tmp_directory


def SAMreader(handle):
    """A simple iterator to yield SAM rows"""

    fields = [
        'QNAME',  # Query NAME of the read or the read pair
        'FLAG',  # Bitwise FLAG (pairing, strand, mate strand, etc.)
        'RNAME',  # Reference sequence NAME
        'POS',  # 1-Based leftmost POSition of clipped alignment
        'MAPQ',  # MAPping Quality (Phred-scaled)
        'CIGAR',  # Extended CIGAR string (operations: MIDNSHP)
        'MRNM',  # Mate Reference NaMe ( if same as RNAME)
        'MPOS',  # 1-Based leftmost Mate POSition
        'ISIZE',  # Inferred Insert SIZE
        'SEQ',  # Query SEQuence on the same strand as the reference
        'QUAL',  # Query QUALity (ASCII-33=Phred base quality)]
    ]

    lines = ifilter(lambda x: not x.startswith('@'), handle)
    for row in csv.reader(lines, delimiter='\t'):
        yield dict(zip(fields, row))


def call_lastz(query_file, ref_file, out_file, outformat='sam'):
    """Calls lastz with the relevant parameters."""

    lastz_path = '/home/will/lastz-distrib/bin/lastz'

    cmd = ' '.join([lastz_path, ref_file, query_file, '--output='+out_file, '--format='+outformat])
    check_call(shlex.split(cmd))


def align_with_lastz(input_seqs, ref_seqs):
    """Aligns set of query sequences with a reference."""

    with tmp_directory() as tmp_dir:
        seq_file = os.path.join(tmp_dir, 'query.fasta')
        ref_file = os.path.join(tmp_dir, 'ref.fasta')
        out_file = os.path.join(tmp_dir, 'res.fasta')

        with open(seq_file, 'w') as handle:
            fasta_writer(handle, input_seqs)

        with open(ref_file, 'w') as handle:
            fasta_writer(handle, ref_seqs)

        call_lastz(seq_file, ref_file, out_file)

        with open(out_file) as handle:
            return list(SAMreader(handle))
