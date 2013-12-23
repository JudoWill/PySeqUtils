""" This tool is for locatating and translating HIV sequences from
 arbitrary locations. It uses multiprocesing to query the LANL
 Sequence Locator tool (http://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html)
 and retrieve the results.

 This tool can either be a library or a command-line tool. If using
 as a command line tool then something like:

  python HIVTransTool.py -o /path/to/output path/to/input

 should do everything you need. To use as a library you should
 import the 'process_seqs' function since it takes care of the
 multiprocesing.

  from HIVTransTool import process_seqs
  for row in process_seqs(input_seqs, threads=5):
    print row

 input_seqs should be an iterable of (name, seq) tuples.

"""
from __future__ import division
__author__ = 'will'
import os.path
from bs4 import BeautifulSoup
from StringIO import StringIO
from GeneralSeqTools import fasta_writer, fasta_reader
from itertools import chain, islice, product, imap
from functools import partial
import glob
import csv
import cookielib
import mechanize
import re
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from copy import deepcopy
from fileinput import FileInput
import argparse
import logging
import time
from httplib import IncompleteRead




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
    prot = None
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


def map_seqs_to_ref(input_seqs, retry=0):
    """Maps a set of (name, seq) pairs to HXB2 using LANL"""

    base_seqs = StringIO()
    fasta_writer(base_seqs, input_seqs)
    base_seqs.seek(0)
    fasta_seqs = base_seqs.read()

    br = build_browser()
    br.open('http://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html')
    logging.debug('Opened Browser to LANL')

    br.select_form(nr=1)
    br.form['SEQ'] = fasta_seqs
    resp = br.submit()
    logging.info('Submitted Seqs to LANL')

    try:
        soup = BeautifulSoup(resp)
    except IncompleteRead:
        if retry > 5:
            raise ValueError
        return map_seqs_to_ref(input_seqs, retry=retry+1)
    rows = []
    count = 0
    for name, seq, table in zip(yield_seq_names(soup), yield_query_seqs(soup), yield_seq_tables(soup)):
        count += 1
        for row in yield_row_vals(table, seq):
            row['Name'] = name
            rows.append(row)

    logging.info('LANL returned %i regions for %i patients' % (len(rows), count))
    return rows


def take(iterable, chunk):
    return list(islice(iterable, chunk))


def yield_chunks(iterable, chunksize):
    """Yields chunks of items from an iterable."""

    chunk = take(iterable, chunksize)
    while chunk:
        yield chunk
        chunk = take(iterable, chunksize)


def process_seqs(input_seqs, threads=5, extract_regions=False, known_names=None):
    """Calls map_seqs_to_ref in a multithreaded way."""

    if extract_regions:
        region_dict = load_region()
    else:
        region_dict = defaultdict(list)

    chunksize = 50
    iter_seqs = iter(input_seqs)
    chunk_iterable = yield_chunks(iter_seqs, chunksize)

    if threads > 1:
        logging.warning('Started ThreadPool with %i workers' % threads)
        ex = ThreadPoolExecutor(max_workers=threads)
        res_iter = chain.from_iterable(ex.map(map_seqs_to_ref, chunk_iterable))
    else:
        logging.warning('Running with NO THREADS!')
        res_iter = chain.from_iterable(imap(map_seqs_to_ref, chunk_iterable))

    name_count = 0
    prev_name = None
    for row in res_iter:
        if row['Name'] != prev_name:
            prev_name = row['Name']
            name_count += 1
        yield row
        for region_row in region_dict[row['RegionName']]:
            nrow = region_linker(deepcopy(row), region_row)
            if nrow:
                yield nrow


def region_linker(row, region_row):

    region_name, start, end, is_aa = region_row

    row['RegionName'] = region_name
    if is_aa:
        out_seq = extract_region(row['QueryAA'],
                                 row['RegionAAStart'],
                                 row['RegionAAStop'],
                                 start, end)
    else:
        out_seq = extract_region(row['QueryNuc'],
                                 row['QueryNucStart'],
                                 row['QueryNucStop'],
                                 start, end)
    if len(out_seq) == 0:
        return None
    if (sum(1 for l in out_seq if l != '-')/len(out_seq)) > 0.5:
        if is_aa:
            row['QueryAA'] = out_seq
        else:
            row['QueryNuc'] = out_seq
        reset_fields = ['QueryNucStart', 'QueryNucStop',
                        'RegionNucStart', 'RegionNucStop',
                        'RegionAAStart', 'RegionAAStop']
        for f in reset_fields:
            row[f] = None
        return row
    return None


def load_region(path='/home/will/PySeqUtils/HIVDBFiles/HXB2RegionNames.csv'):

    region_dict = defaultdict(list)
    with open(path) as handle:
        for row in csv.DictReader(handle):
            region_dict[row['SourceName']].append((
                row['RegionName'],
                int(row['RegionStart'])-1,  # since the file is ONE-BASED
                int(row['RegionEnd'])-1,
                row['IsAA'] == 'True'
            ))
    return region_dict


def extract_region(found_seq, found_start, found_stop, region_start, region_stop, fillval='-'):
    """Extracts strings from the region. ZERO-BASED INDEXING!!"""

    if found_start > region_start:
        found_seq = fillval*abs(region_start - found_start) + found_seq
        found_start = region_start
    if found_stop < region_stop:
        found_seq = found_seq + fillval*abs(found_stop - region_stop)

    start = (region_start - found_start)
    end = start + (region_stop - region_start)
    return found_seq[start:end]


def main(input_files, out_csv_file,
         threads=5, extract_regions=True, known_names=None):

    seq_iter = fasta_reader(FileInput(input_files))
    found_names = set()
    need_header = True
    if os.path.exists(out_csv_file):
        need_header = False
        with open(out_csv_file) as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                found_names.add(row['Name'])

    csv_handle = open(out_csv_file, 'a')
    fields = ['Name', 'RegionName', 'QueryNucStart', 'QueryNucStop', 'QueryNuc',
              'RegionNucStart', 'RegionNucStop', 'RegionAAStart', 'RegionAAStop', 'QueryAA']
    csv_writer = csv.DictWriter(csv_handle, fields, delimiter='\t')
    if need_header:
        csv_writer.writeheader()

    wanted_seqs = ((name, seq) for name, seq in seq_iter if name not in found_names)
    result_iter = process_seqs(wanted_seqs,
                               threads=threads,
                               extract_regions=extract_regions,
                               known_names=known_names)
    for chunk in yield_chunks(result_iter, 1000):
        print chunk[0]
        csv_writer.writerows(chunk)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A tool for extracting reading frames from HIV sequences')

    parser.add_argument('infiles', type=str, help='Input files in Fasta format')
    parser.add_argument('-o', type=str, required=True, help='Output template.')
    parser.add_argument('-R', action='store_true', default=True, help='Extract internal regions like V3?')
    parser.add_argument('-t', type=int, default=5, help='Number of threads to use when querying LANL. DONT BE A DICK!')
    parser.add_argument('-q', action='store_true', default=False, help='Be Quiet!')

    args = parser.parse_args()

    fmt = '%(threadName)s: %(levelname)s %(asctime)s | %(message)s'
    if args.q:
        logging.basicConfig(level=logging.WARNING, format=fmt)
    else:
        logging.basicConfig(level=logging.INFO, format=fmt)

    infiles = glob.glob(args.infiles)
    out_csv = args.o
    nthreads = args.t

    logging.info('Calculating the number of seqs in %s' % ','.join(infiles))
    known_names = sum(line.startswith('>') for line in FileInput(infiles))
    logging.warning('Found %i sequences to process' % known_names)
    logging.warning('Starting!')
    main(infiles, out_csv,
         threads=nthreads, extract_regions=args.R, known_names=known_names)
    logging.warning('Finished!')