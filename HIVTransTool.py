__author__ = 'will'

from StringIO import StringIO
from GeneralSeqTools import fasta_writer
from QuasiTools import SAMreader
from subprocess import check_output
from tempfile import NamedTemporaryFile as NTF
import shlex


def map_seqs_to_ref(input_seqs, ref_path):
    """Maps a set of (name, seq) pairs to the reference. Returns the start-position and the mapped sequence."""

    with NTF('w') as in_seqs:

        fasta_writer(in_seqs, input_seqs)
        in_seqs.seek(0)

        cmd = 'lastz %s --ambiguous=iupac --format=sam'
        ncmd = cmd % ref_path

        out = check_output(shlex.split(ncmd), stdin=in_seqs)
        for row in SAMreader(StringIO(out)):
            yield row['QNAME'], row['POS'], row['SEQ']
