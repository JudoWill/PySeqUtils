__author__ = 'will'
from itertools import groupby, islice, imap
import subprocess # import check_output
from StringIO import StringIO
from tempfile import NamedTemporaryFile as NTF
from concurrent.futures import ThreadPoolExecutor
import shlex
import os
from pandas import DataFrame
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import IUPAC


def write_nexus_alignment(input_seqs, handle, is_aa=True):

    alphabet = IUPAC.protein if is_aa else IUPAC.unambiguous_dna
    seqs = []
    tmp_handle = StringIO()
    for name, seq in input_seqs:
        nseq = ''.join(seq).replace('O', '-')
        bseq = SeqRecord(Seq(nseq, alphabet=alphabet), id=name)
        seqs.append(bseq)
    SeqIO.write(seqs, tmp_handle, 'nexus')
    tmp_handle.seek(0)
    strdata = tmp_handle.read().replace("'", '')
    handle.write(strdata)


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

    return list(fasta_reader(StringIO(out)))


def seq_map_to_ref(seq_align, ref_align):

    out_seq = []
    for s, r in zip(seq_align, ref_align):
        if r != '-':
            out_seq.append(s)
    return ''.join(out_seq)


def seq_align_to_ref(input_seqs, ref_seq, max_workers=None):
    """Aligns all sequences to a reference.
    """

    check_seqs = [[(name, seq), ('__ref__', ref_seq)] for name, seq in input_seqs]
    if max_workers > 1:
        executor = ThreadPoolExecutor(max_workers=max_workers)
        res = executor.map(call_muscle, check_seqs)
    else:
        res = imap(call_muscle, check_seqs)

    for alignment in res:
        name = alignment[0][0]
        ref_align = seq_map_to_ref(alignment[0][1], alignment[1][1])
        yield name, ref_align


def convert_seqs_to_dataframe(input_seqs):

    items = {}
    for name, seq in input_seqs:
        items[name] = list(seq)

    return DataFrame(items).T


def convert_seqDF_to_list(frame):

    print frame

    out_seqs = []
    for key, row in frame.iterrows():
        out_seqs.append((key, ''.join(row)))

    return out_seqs


def extract_region_from_alignment(input_aln, start, stop):
    """Extracts a region from an alignment
    """

    for name, seq in input_aln:
        yield name, seq[start:stop]

