from __future__ import division
__author__ = 'will'
import os
from bs4 import BeautifulSoup
from StringIO import StringIO
from GeneralSeqTools import fasta_writer
from itertools import chain, islice
import csv
import cookielib
import mechanize
import re
from concurrent.futures import ThreadPoolExecutor


def build_browser():
    """Builds a browser which can 'fool' the LANL database."""

    br = mechanize.Browser()

    # Cookie Jar
    cj = cookielib.LWPCookieJar()
    br.set_cookiejar(cj)

    # Browser options
    br.set_handle_equiv(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)

    # Follows refresh 0 but not hangs on refresh > 0
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)

    # User-Agent (this is cheating, ok?)
    br.addheaders = [('User-agent', 'Mozilla/5.0 (Windows; U; Windows NT 6.0; en-US; rv:1.9.0.6')]

    return br


def get_nums(instr):
    """Returns a list of integers that occur in text."""
    return map(int, re.findall('\d+', instr))


def is_seq_table(tab):
    """Determines whether a BS4 object is a table which has sequence results."""

    wanted_str = 'Table of genomic regions touched by query sequence. (Protein translation of query shown in blue.)'
    check_str = ''.join(l for l in wanted_str if not l.isspace())
    frow = tab.tr.td
    if frow is None:
        return False
    nstr = ''.join(l for l in frow.text if not l.isspace())
    if nstr == check_str:
        return True
    else:
        return False


def yield_query_seqs(soup):
    """Returns the query sequence associate with a particular query."""

    for table in soup.findAll('table'):
        if table.form and table.form.input.attrs['name'] == 'PasteQuery':
            yield table.form.input.attrs['value']


def yield_seq_tables(soup):
    """Yields a table which has translated sequence information."""

    for table in soup.findAll('table'):
        if is_seq_table(table):
            yield table


def yield_seq_names(soup):
    """Yields the query name from a LANL output."""

    for hd in soup.findAll('h3'):
        yield hd.text.replace('Sequence ', '')


def yield_row_vals(table, nuc_seq):
    """Processes a 'sequence table' from LANL and returns the contained regions."""

    is_seq_row = False
    for row in table.findAll('tr'):
        cols = list(row.findAll('td'))
        #print cols
        if len(cols) == 0:
            is_seq_row = True
        if is_seq_row and (len(cols) == 5):
            prot = cols[0].text.strip()
            hxb2_nuc_rel = get_nums(cols[1].text.strip())
            query_nuc_rel = get_nums(cols[2].text.strip())
            hxb2_aa_rel = get_nums(cols[4].text.strip())
            if 'LTR' in prot:
                odict = {
                    'RegionName': prot,
                    'QueryNucStart': query_nuc_rel[0],
                    'QueryNucStop': query_nuc_rel[1],
                    'QueryNuc': nuc_seq[query_nuc_rel[0]-1:query_nuc_rel[1]-1],
                    'RegionNucStart': hxb2_nuc_rel[0],
                    'RegionNucStop': hxb2_nuc_rel[1],
                    'RegionAAStart': None,  #hxb2_aa_rel[0],
                    'RegionAAStop': None,  #hxb2_aa_rel[1],
                    'QueryAA': None
                }
                yield odict

        elif is_seq_row and (len(cols) == 1) and ('Notice' not in cols[0].text):
            seq = ''.join(l for l in cols[0].text if not l.isspace())
            odict = {
                'RegionName': prot,
                'QueryNucStart': query_nuc_rel[0],
                'QueryNucStop': query_nuc_rel[1],
                'QueryNuc': nuc_seq[query_nuc_rel[0]-1:query_nuc_rel[1]-1],
                'RegionNucStart': hxb2_nuc_rel[0],
                'RegionNucStop': hxb2_nuc_rel[1],
                'RegionAAStart': hxb2_aa_rel[0],
                'RegionAAStop': hxb2_aa_rel[1],
                'QueryAA': seq
            }
            yield odict


def map_seqs_to_ref(input_seqs):
    """Maps a set of (name, seq) pairs to HXB2 using LANL"""

    base_seqs = StringIO()
    fasta_writer(base_seqs, input_seqs)
    base_seqs.seek(0)
    fasta_seqs = base_seqs.read()

    br = build_browser()
    br.open('http://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html')

    br.select_form(nr=1)
    br.form['SEQ'] = fasta_seqs
    resp = br.submit()

    soup = BeautifulSoup(resp)
    rows = []
    for name, seq, table in zip(yield_seq_names(soup), yield_query_seqs(soup), yield_seq_tables(soup)):
        for row in yield_row_vals(table, seq):
            row['Name'] = name
            rows.append(row)
    return rows


def take(iterable, chunk):
    return list(islice(iterable, chunk))


def yield_chunks(iterable, chunksize):
    """Yields chunks of items from an iterable."""

    chunk = take(iterable, chunksize)
    while chunk:
        yield chunk
        chunk = take(iterable, chunksize)


def process_seqs(input_seqs, threads=5):
    """Calls map_seqs_to_ref in a multithreaded way."""

    chunksize = 10
    iter_seqs = iter(input_seqs)

    if threads > 1:
        ex = ThreadPoolExecutor(max_workers=threads)
        process = ex.map
    else:
        process = map

    chunk_iterable = yield_chunks(iter_seqs, chunksize)
    res_iter = chain.from_iterable(process(map_seqs_to_ref, chunk_iterable))
    for row in res_iter:
        yield row

