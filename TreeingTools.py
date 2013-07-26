__author__ = 'will'

import os
import dendropy
import contextlib
from subprocess import check_output
import shlex
from tempfile import mkdtemp
from tempfile import NamedTemporaryFile as NTF
import shutil
from GeneralSeqTools import write_nexus_alignment
import csv
from StringIO import StringIO
from dendropy.treecalc import PatristicDistanceMatrix
from itertools import combinations
from collections import defaultdict
from scipy.stats import ttest_ind
from random import shuffle
import numpy as np
import logging
from types import ListType, StringType
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymongo import MongoClient

from celery import Celery, subtask
from celery import group

celery = Celery()
celery.config_from_object('celeryconfig')



@celery.task(name='TreeingTools.write_results_to_mongo',
             queue='writingqueue')
def write_results_to_mongo(result, **kwargs):
    #, result_type=None, database=None, extra_fields=None

    if type(result) == ListType:
        result_list = list(result)
    else:
        result_list = [result]

    client = MongoClient()
    db = client['new_test_db'].results
    for res in result_list:
        if kwargs['extra_fields']:
            res.update(kwargs['extra_fields'])
        res.pop(None, None)
        res['ResultType'] = kwargs['result_type']
        db.insert(res)


@celery.task(name='TreeingTools.process_region',
             queue='base')
def process_region(input_seqs, trop_dict, basename, database=None, mrbayes_args={}, extra_fields={}):

    try:
        if os.path.exists(basename + '.tree'):
            return
        handle = open(basename + '.tree', 'w')
        logging.info('Making Tree ' + str(extra_fields))
        contree, treeset = make_mrbayes_trees(input_seqs, **mrbayes_args)
        contree.write_to_stream(handle, 'nexus')
        treeset.write_to_path(basename + '.treeset', 'nexus')
    except IOError:
        return
    except OSError:
        return

    bats_write_subtask = subtask('TreeingTools.write_results_to_mongo', (), {
        'result_type': 'BATS',
        'extra_fields': extra_fields,
        'database': database
    })

    logging.info('Starting BATS ' + str(extra_fields))
    run_bats.apply_async(args=(basename + '.treeset', trop_dict),
                         kwargs = {'nreps': 50}, link = bats_write_subtask)

    benj_write_subtask = subtask('TreeingTools.write_results_to_mongo', (), {
        'result_type': 'Benj',
        'extra_fields': extra_fields,
        'database': database
    })

    dmat = get_pairwise_distances(contree)
    logging.info('Starting Dist Pvals ' + str(extra_fields))
    check_distance_pvals.apply_async(args=(dmat, trop_dict),
                                     kwargs = {'nreps': 500}, link = benj_write_subtask)



@contextlib.contextmanager
def tmp_directory(rm_dir=True, *args, **kwargs):
    """A context manager which changes the working directory to the given
    path, and then changes it back to its previous value on exit.

    """
    path = mkdtemp(*args, **kwargs)
    try:
        yield path + '/'
    finally:
        if rm_dir:
            shutil.rmtree(path)

@contextlib.contextmanager
def push_dir(path):

    cur_path = os.path.abspath(os.curdir)
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cur_path)


def clean_sequences(input_seqs, alphabet=None, is_aa=True):

    if (alphabet == generic_protein) or is_aa:
        allowed = set(IUPAC.IUPACProtein.letters)
    elif (alphabet == generic_dna) or not is_aa:
        allowed = set(IUPAC.IUPACUnambiguousDNA.letters)
    else:
        raise KeyError('Unknown alphabet!')

    for name, seq in input_seqs:
        nseq = ''.join(l if l.upper() in allowed else '-' for l in seq)
        yield name, nseq


@celery.task(name='TreeingTools.make_mrbayes_trees',
             queue='long-running')
def make_mrbayes_trees(input_seqs, mrbayes_kwargs=None, is_aa=True):
    """Takes an ALIGNED set of sequences and generates a phylogenetic tree using MrBayes.
    """

    cleaned_seqs = clean_sequences(input_seqs, is_aa=is_aa)
    with tmp_directory() as tmpdir:
        align_file = tmpdir + 'seqalign.nxs'
        mrbayes_cmd_file = tmpdir + 'analysis.nxs'
        multi_prob = tmpdir + 'seqalign.nxs.trprobs'
        cons_file = tmpdir + 'seqalign.nxs.con.tre'

        with open(align_file, 'w') as handle:
            write_nexus_alignment(cleaned_seqs, handle, is_aa=is_aa)

        with open(mrbayes_cmd_file, 'w') as handle:
            if mrbayes_kwargs is None:
                mrbayes_kwargs = {}
            txt = generate_mrbayes_nexus(align_file, align_file, is_aa=is_aa, **mrbayes_kwargs)
            handle.write(txt)

        cmd = 'mb ' + mrbayes_cmd_file
        print cmd
        check_output(shlex.split(cmd))

        with open(multi_prob) as handle:
            trees = dendropy.TreeList.get_from_stream(handle, schema='nexus')
        with open(cons_file) as handle:
            con_tree = dendropy.Tree.get_from_stream(handle, schema='nexus')

    return con_tree, trees


def generate_mrbayes_nexus(alignment_path, output_path,
                           nchains=1, ngen=5000, samplefreq=1000,
                           is_aa=True):
    """Generates the NEXUS command to control MrBayes in the form that I usually use. This will likely be expanded as I
     have to include more issues.
    """

    cmd = """begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute %(align)s;
   %(model)s;
   mcmc nchains = %(nchains)i ngen = %(ngen)i samplefreq=%(samplefreq)i file=%(out)s;
   sump;
   sumt;
end;"""

    tdict = {
        'align': alignment_path,
        'out': output_path,
        'model': 'prset aamodelpr = mixed' if is_aa else 'lset nst=6 rates=invgamma',
        'nchains': nchains,
        'ngen': ngen,
        'samplefreq': samplefreq
    }

    return cmd % tdict


def bats_format_nexus(treeset, outhandle, trop_dict):
    """Writes the treeset in the format required by BatS.
    """

    def process_tree_line(line):

        parts = line.strip().split()
        return 'tree %s [&R] %s\n' % (parts[1], parts[-1])

    def print_items(taxon):
        return conv_dict.get(str(taxon), None)

    leafs = set()
    for tree in treeset:
        for leaf in tree.leaf_iter():
            leafs.add(str(leaf.taxon))
    leafs = sorted(leafs, key=lambda x: len(x), reverse=True)

    outhandle.write('#NEXUS\n\n\n')
    outhandle.write('begin states;\n')
    conv_dict = {}
    for num, leaf in enumerate(leafs, 1):
        outhandle.write('%i %s\n' % (num, trop_dict[leaf]))
        conv_dict[leaf] = str(num)

    outhandle.write('End;\n\n')

    outhandle.write('begin trees;\n')
    for num, tree in enumerate(treeset, 1):
        tstr = tree.as_newick_string(reverse_translate=print_items, suppress_edge_lengths=True)
        outhandle.write('tree tree_%i [&R] %s;\n' % (num, tstr))
    outhandle.write('end;\n')

@celery.task(name='TreeingTools.run_bats',
             queue='long-running')
def run_bats(treeset, trop_dict, nreps=5000):
    """Runs the BatS analysis on a treeset.
    """

    if type(treeset) is StringType:
        treeset = dendropy.TreeList.get_from_path(treeset, 'nexus')
    nstates = len(set(trop_dict.values()))
    with NTF() as handle:
        bats_format_nexus(treeset, handle, trop_dict)
        handle.flush()
        os.fsync(handle)
        cmd = 'java -jar /home/will/BaTS_beta_build2.jar single %s %i %i'
        logging.info('Running BATS')
        out = check_output(shlex.split(cmd % (handle.name, nreps, nstates)))
        logging.info('Sucessful BATS')

    handle = StringIO(out)
    headers = []
    for line in handle:
        if line.startswith('Stat'):
            headers = line.strip().split('\t')
            break
    return list(csv.DictReader(handle, fieldnames=headers, delimiter='\t'))[:-2]

@celery.task(name='TreeingTools.get_pairwise_distances')
def get_pairwise_distances(tree):
    """Returns a dict of pairwise distances."""

    taxons = tree.taxon_set
    pdm = PatristicDistanceMatrix(tree)
    pdm.calc()
    dmat = {}
    for p1, p2 in combinations(taxons, 2):
        d = pdm(p1, p2)
        dmat[(p1.label, p2.label)] = d
        dmat[(p2.label, p1.label)] = d

    return dmat

@celery.task(name='TreeingTools.check_distance_pvals')
def check_distance_pvals(dist_dict, group_dict, group_frac=0.5, nreps=500):

    groups = sorted(set(group_dict.values()))
    assert len(groups) == 2

    group_vals = defaultdict(list)

    for (key1, key2), dist in dist_dict.items():
        if group_dict[key1] == group_dict[key2]:
            group_vals[group_dict[key1]].append(dist)

    _, raw_pval = ttest_ind(*group_vals.values())

    nitems = int(group_frac*min(map(len, group_vals.values())))
    cor_vals = []
    for _ in range(nreps):
        [shuffle(items) for items in group_vals.values()]
        _, pval = ttest_ind(*[items[:nitems] for items in group_vals.values()])
        cor_vals.append(pval)

    odict = {
        'RawPval': raw_pval,
        'AdjPval': np.mean(cor_vals),
        'Group1Name': groups[0],
        'Group2Name': groups[1],
        'Group1Mean': np.mean(group_vals[groups[0]]),
        'Group2Mean': np.mean(group_vals[groups[1]]),
        'Group1Std': np.std(group_vals[groups[0]]),
        'Group2Std': np.std(group_vals[groups[1]])
    }

    return odict
