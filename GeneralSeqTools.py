__author__ = 'will'
from itertools import groupby, islice
import subprocess # import check_output
from StringIO import StringIO
from tempfile import NamedTemporaryFile as NTF
import shlex
import os


def fasta_reader(handle):
    """Reads a fasta formatted handle and yields (name, seq) pairs
    """

    name = None
    for key, lines in groupby(handle, key=lambda x: x.startswith('>')):
        if key:
            name = lines.next()[1:].strip()
        else:
            seq = ''.join(line.strip() for line in lines)
            yield name, seq


def take(num, iterable):
    return list(islice(iterable, num))


def fasta_writer(handle, seqs):
    """Writes (name, seq) pairs in fasta format
    """

    width = 80
    for name, seq in seqs:
        handle.write('>%s\n' % name)
        siter = iter(seq)
        block = take(width, siter)
        while block:
            handle.write('%s\n' % ''.join(block))
            block = take(width, siter)


def call_muscle(input_seqs):
    """Calls muscle and aligns the input sequences.
    """

    with NTF() as seq_handle:
        fasta_writer(seq_handle, input_seqs)
        seq_handle.flush()
        os.fsync(seq_handle.fileno())
        cmd = 'muscle -in %s -quiet' % seq_handle.name
        out = subprocess.check_output(shlex.split(cmd))

    return fasta_reader(StringIO(out))








