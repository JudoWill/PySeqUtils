from __future__ import division
__author__ = 'will'

import os
import dendropy
import contextlib
from subprocess import check_output, check_call
import shlex
from tempfile import mkdtemp
from tempfile import NamedTemporaryFile as NTF
import shutil
from GeneralSeqTools import write_nexus_alignment, fasta_writer
import csv
from StringIO import StringIO
from dendropy.treecalc import PatristicDistanceMatrix
from itertools import combinations, product
from collections import defaultdict, Counter
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
from scipy.stats import gaussian_kde

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


def clean_sequences(input_seqs, alphabet=None, is_aa=True, dump_if_all_gaps=True):

    if (alphabet == generic_protein) or is_aa:
        allowed = set(IUPAC.IUPACProtein.letters)
    elif (alphabet == generic_dna) or not is_aa:
        allowed = set(IUPAC.IUPACUnambiguousDNA.letters)
    else:
        raise KeyError('Unknown alphabet!')

    for name, seq in input_seqs:
        nseq = ''.join(l if l.upper() in allowed else '-' for l in seq.upper())
        if all(l == '-' for l in nseq) and dump_if_all_gaps:
            continue

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
                           nchains=2, ngen=50000, samplefreq=10000,
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

    assert len(group_vals) == 2
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


def write_phylip_seqs(seqs, outhandle, alphabet=generic_protein):
    """Writes sequences to a Phylip formatted file. Due to the nature of the
    format it standarizes the names and returns a dict mapping the orignal name
    to the new name.

    seqs -- A sequence of (name, seq) tuples.
    outhandle -- The handle to write sequences into.
    alphabet = generic_protein -- This is to define the input type of the
      sequences.

    returns:
    trans_names -- A dict() mapping orignal names to the standardized names.
    """

    trans_names = {}
    out_seqs = []
    for num, (name, seq) in enumerate(clean_sequences(seqs, alphabet=alphabet)):
        new_name = 'Seq-%i' % num
        trans_names[name] = new_name
        out_seqs.append(SeqRecord(Seq(seq,
                                      alphabet=alphabet),
                                  id=str(new_name)))
    SeqIO.write(out_seqs, outhandle, 'phylip')
    return trans_names


def make_phylip_seq_dist_mat(inseqs, alphabet, tmp_path=None, rm_dir=True):
    """Takes the input sequences and uses the PHYLIP tool to find their
    pairwise distance. This is then used to do the neighbor joining in
    later analysis.

    inseqs -- A sequence of (name, sequence) tuples.
    alphabet -- A Bio.Alphabet object indicating the the sequence type.
                Currently this only accepts generic_protein and
                generic_dna
    tmp_path=None -- A path to create the temporary directories.
    rm_dir=True -- Whether to remove the temporary directory after completion.

    returns:
    trans_names -- A dict() mapping the new-names to the standardized names.
    dist_mat -- A str of the data returned from the phylip protdist
                function. This is not programatically useful but can
                be written to a file for later analysis.

    """

    commands = "2\ny\n"

    with tmp_directory(dir=tmp_path, rm_dir=rm_dir) as r_path:
        with push_dir(r_path):
            with open('cmd_file', 'w') as cmd_handle:
                cmd_handle.write(commands)
            cmd_handle = open('cmd_file')
            with open('infile', 'w') as handle:
                trans_names = write_phylip_seqs(inseqs, handle, alphabet = alphabet)
            if alphabet == generic_protein:
                cmd = 'phylip protdist'
            elif alphabet == generic_dna:
                cmd = 'phylip dnadist'
            else:
                raise KeyError('Unknown alphabet.')
            cmd = shlex.split(cmd)
            check_call(cmd, stdin=cmd_handle)

            with open('outfile') as handle:
                dist_data = handle.read()
                return trans_names, dist_data


def process_phylip_dist_mat(dist_data, trans_names):
    """Converts the output of make_phylip_seq_dist_mat into a format
    appropriate to use in downstream analysis.

    dist_data       A string of the dist_data as returned from the
                    phylip protdist and dnadist programs.
    trans_names     A dict mapping the translated names used in the
                    phylip alignment format with the original names.

    Returns:
    dmat -- A dict() keyed by all pairs of sequences with the genetic distance
            as the value. This can be used in the check_distance_pvals
            function.
    """

    rev_dict = dict((val, key) for key, val in trans_names.items())
    handle = StringIO(dist_data)
    handle.next()
    dist_groups = []

    for line in handle:
        if line.startswith('Seq-'):
            dist_groups.append('')
        dist_groups[-1] += line

    omat = {}
    tmpl = 'Seq-%i'
    for seq_num, group in enumerate(dist_groups):
        parts = group.split()[1:]
        for onum, val in enumerate(parts):
            nkey = (rev_dict[tmpl % seq_num], rev_dict[tmpl % onum])
            nval = float(val)
            if nval >= 0:
                omat[nkey] = nval

    return omat


def make_phylip_tree(dist_data, tmp_path=None, rm_dir=True):
    """Takes the dist_data returned by make_phylip_seq_dist_mat the
    phylip neighbor tool to generate a tree.

    dist_data -- The distance data returned by make_phylip_seq_dist_mat.
    tmp_path=None -- A path to create the temporary directories.
    rm_dir=True -- Whether to remove the temporary directory after completion.

    returns:
    tree -- A dendropy tree object of the resulting tree.
    """

    commands = "2\n3\ny\n"
    with tmp_directory(dir=tmp_path, rm_dir=rm_dir) as r_path:
        with push_dir(r_path):
            with open('cmd_file', 'w') as cmd_handle:
                cmd_handle.write(commands)
            cmd_handle = open('cmd_file')
            with open('infile', 'w') as handle:
                handle.write(dist_data)
            cmd = shlex.split('phylip neighbor')
            check_call(cmd, stdin = cmd_handle)
            with open('outtree') as handle:
                tree = dendropy.Tree.get_from_stream(handle, schema='newick')
                return tree


def phylip_tree(seqs, alphabet=generic_protein, tmp_path=None, rm_dir=True):
    """Uses the Phylip program to generate a tree from the input sequences.

    seqs -- A sequence of (name, sequence) tuples.
    alphabet -- A Bio.Alphabet object indicating the the sequence type.
                Currently this only accepts generic_protein and
                generic_dna
    tmp_path=None -- A path to create the temporary directories.
    rm_dir=True -- Whether to remove the temporary directory after completion.

    returns:
    tree -- A dendropy tree object of the resulting tree.
    dist_data -- The a dict() of phylogentic distances between all pairs.
    """

    trans_names, dist_data = make_phylip_seq_dist_mat(seqs, alphabet,
                                                      tmp_path=tmp_path,
                                                      rm_dir=rm_dir)
    out_tree = make_phylip_tree(dist_data,
                                tmp_path=tmp_path,
                                rm_dir=rm_dir)
    for orig_name, new_name in trans_names.items():
        node = out_tree.find_node_with_taxon_label(new_name)
        if node:
            node.taxon.label = orig_name
    return out_tree, process_phylip_dist_mat(dist_data, trans_names)


def run_FastTree(seqs, alphabet=generic_protein, tmp_path=None, uniq_seqs=False):

    if uniq_seqs:

        trans_names = defaultdict(list)
        norm_seq_names = {}
        for num, (name, seq) in enumerate(seqs):
            trans_names[seq].append(name)
            new_name = 'Seq-%i' % num
            norm_seq_names[seq] = new_name

        uni_seqs = []
        name_defs = {}
        for seq, new_name in norm_seq_names.items():
            uni_seqs.append((new_name, seq))
            name_defs[new_name] = trans_names[seq]

        out_tree = run_FastTree(uni_seqs, alphabet=alphabet, tmp_path=tmp_path)

        tax_set = out_tree.taxon_set
        for old_name, new_names in name_defs.items():
            node = out_tree.find_node_with_taxon_label(old_name)
            if node:
                names = iter(new_names)
                node.taxon.label = names.next()
                parent = node.parent_node
                edge_dist = node.edge.length
                for name in names:
                    parent.new_child(taxon=tax_set.new_taxon(label=name),
                                     edge_length=edge_dist)
        return out_tree

    else:
        cmd = '/home/will/PySeqUtils/FastTree %(alpha)s -quiet %(path)s'

        with NTF(dir=tmp_path, suffix='.fasta') as handle:
            fasta_writer(handle, seqs)
            handle.flush()
            os.fsync(handle)
            tdict = {
                'alpha': '-nt' if alphabet == generic_dna else '',
                'path': handle.name
            }
            cmd_list = shlex.split(cmd % tdict)
            tree_str = check_output(cmd_list)
        return dendropy.Tree(stream=StringIO(tree_str), schema='newick')


def fast_tree(seqs, alphabet=generic_protein, tmp_path=None, rm_dir=True):
    """A dropin replacement for phylip that is significantly faster for larger
     sequence alignments.

    seqs -- A sequence of (name, sequence) tuples.
    alphabet -- A Bio.Alphabet object indicating the the sequence type.
                Currently this only accepts generic_protein and
                generic_dna
    tmp_path=None -- A path to create the temporary directories.
    rm_dir=True -- Whether to remove the temporary directory after completion. **Does NOTHING**

    returns:
    tree -- A dendropy tree object of the resulting tree.
    dist_data -- The a dict() of phylogentic distances between all pairs.
    """

    out_tree = run_FastTree(seqs, alphabet=alphabet, tmp_path=tmp_path)
    mat_calc = dendropy.treecalc.PatristicDistanceMatrix(out_tree)
    dmat = {}
    for p1, p2 in combinations(out_tree.taxon_set, 2):
        dist = mat_calc(p1, p2)
        dmat[(p1.label, p2.label)] = dist
        dmat[(p2.label, p1.label)] = dist

    return out_tree, dmat


def phylip_tree_collapse_unique(seqs, alphabet=generic_protein, tmp_path=None, rm_dir=True, use_fast=True):
    """Uses the Phylip program to generate a tree from the input sequences.
    However, this one intelligently deals with repeated sequences.

    seqs -- A sequence of (name, sequence) tuples.
    alphabet -- A Bio.Alphabet object indicating the the sequence type.
                Currently this only accepts generic_protein and
                generic_dna
    tmp_path=None -- A path to create the temporary directories.
    rm_dir=True -- Whether to remove the temporary directory after completion.

    returns:
    tree -- A dendropy tree object of the resulting tree.
    dist_data -- The phylogentic distance between all pairs.
    """

    name_defs = defaultdict(set)
    seq_names = {}

    uni_nseqs = []
    for name, seq in seqs:
        if seq not in seq_names:
            new_name = 'UniSeq-%i' % len(seq_names)
            uni_nseqs.append((new_name, seq))
            seq_names[seq] = new_name
            name_defs[new_name].add(name)
        else:
            name_defs[seq_names[seq]].add(name)

    if len(uni_nseqs) < 4:
        raise ValueError('Too few unique sequences')

    if use_fast:
        out_tree, dist_mat = fast_tree(uni_nseqs, alphabet=alphabet,
                                       tmp_path=tmp_path, rm_dir=rm_dir)
    else:
        out_tree, dist_mat = phylip_tree(uni_nseqs, alphabet=alphabet,
                                         tmp_path=tmp_path, rm_dir=rm_dir)

    tax_set = out_tree.taxon_set
    for old_name, new_names in name_defs.items():
        node = out_tree.find_node_with_taxon_label(old_name)
        if node:
            names = iter(new_names)
            node.taxon.label = names.next()
            parent = node.parent_node
            edge_dist = node.edge.length
            for name in names:
                parent.new_child(taxon=tax_set.new_taxon(label=name),
                                 edge_length=edge_dist)

    new_dmat = {}
    for (n1, n2), dist in dist_mat.items():
        for new_1, new_2 in product(name_defs[n1], name_defs[n2]):
            new_dmat[(new_1, new_2)] = dist

    return out_tree, new_dmat


def calculate_ai(tree, pheno_dict):
    """Calculates the Association Index for a tree and a set of phenotypes.
    http://jvi.asm.org/content/75/23/11686.full

    Starting from the root of the tree, the composition of sequences in
    each successive bifurcating node was calculated. An association value, d,
    for the tree was calculated by summation of values individually
    calculated from each node, according to the formula:

    d = (1-f)/2**(n-1),

    n = # sequences below this node
    f = the frequency of most common sample type.

    tree -- A dendropy tree.
    pheno_dict -- A dict() of phenotypes.

    returns:
    AI
    """

    ai = 0.0
    for node in tree.nodes():
        if ~node.is_leaf():
            leafs = node.leaf_nodes()
            counts = Counter(pheno_dict[leaf.taxon.label] for leaf in leafs)
            best_count = counts.most_common(n=1)[0][1]
            f = best_count/len(leafs)

            #have to do it in log-space due to Overflow errors with large trees
            ai += 2**(np.log2(1-f)-(len(leafs)-1))

    return ai


def evaluate_association_index(tree, pheno_dict, nreps=100):

    true_ai = calculate_ai(tree, pheno_dict)

    rand_vals = []
    new_tree = dendropy.Tree(tree)
    for _ in range(nreps):
        new_tree.randomly_assign_taxa(new_tree.taxon_set)
        rand_vals.append(calculate_ai(new_tree, pheno_dict))
    kde = gaussian_kde(rand_vals)
    pval = kde.integrate_box_1d(true_ai, np.inf)
    rand_mean = np.mean(rand_vals)
    return true_ai, pval, rand_mean