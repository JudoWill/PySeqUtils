__author__ = 'will'
from itertools import groupby, islice, imap
import subprocess  # import check_output
from StringIO import StringIO
from tempfile import NamedTemporaryFile as NTF
from concurrent.futures import ThreadPoolExecutor
import shlex
import os
import re
import csv
import time
from pandas import DataFrame
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import IUPAC
from mechanize import Browser, HTTPError
from collections import defaultdict
from functools import partial
from bs4 import BeautifulSoup


def yield_chunks(iterable, chunksize):

    chunk = take(chunksize, iterable)
    while chunk:
        yield chunk
        chunk = take(chunksize, iterable)


def WebPSSM_V3_series(input_V3_seqs, threads=20, matrix='x4r5'):

    if matrix not in {'x4r5', 'sinsi', 'subC'}:
        raise KeyError('%s matrix not understood.' % matrix)

    if threads:
        executor = ThreadPoolExecutor(max_workers=threads)
        processor = executor.map
    else:
        processor = imap

    pssm_func = partial(WebPSSM_V3_fasta, matrix=matrix)
    chunk_iter = yield_chunks(iter(input_V3_seqs), 30)
    fields = ['name', 'score',
              'pred', ' x4.pct',
              'r5.pct', 'geno',
              'pos.chg', 'net.chg',
              'percentile']
    for res in processor(pssm_func, chunk_iter):
        for row in res:
            yield dict(zip(fields, row))


def WebPSSM_V3_fasta(V3_tuple, matrix='x4r5'):

    if matrix not in {'x4r5', 'sinsi', 'subC'}:
        raise KeyError('%s matrix not understood.' % matrix)

    V3_fasta = ''
    for name, seq in V3_tuple:
        V3_fasta += '>%s\n%s\n' % (name, seq)
    baseurl = 'http://indra.mullins.microbiol.washington.edu/webpssm/'
    basejob = 'http://indra.mullins.microbiol.washington.edu/webpssm/data/results/%s.pssm.txt'
    br = Browser()
    br.open(baseurl)
    br.select_form(nr=0)

    br['seqs'] = V3_fasta
    br['matrix'] = ['subC']
    br['datatype'] = ['aa']
    resp = br.submit()

    soup = BeautifulSoup(resp.read())
    try:
        jobid = re.findall('jobid=(\d+)', soup.meta.attrs['content'])[0]
    except IndexError:
        return []
    joburl = basejob % jobid

    out = None
    attempts = 0
    waittime = 10
    while (out is None) and (attempts < 10):
        try:
            attempts +=1
            nresp = br.open(joburl)
            out = nresp.read()
        except HTTPError:
            time.sleep(waittime)

    if out is None:
        return []
    return list(csv.reader(StringIO(out), delimiter='\t'))[2:]


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


def fasta_writer(handle, seqs, single_line=False):
    """Writes (name, seq) pairs in fasta format
    """

    width = 80
    for name, seq in seqs:
        handle.write('>%s\n' % name)
        siter = iter(seq)
        if single_line:
            width = len(seq) + 1
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
        adict = dict(alignment)
        name = [key for key in adict.keys() if key != '__ref__'][0]
        ref_align = seq_map_to_ref(adict[name], adict['__ref__'])
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

