__author__ = 'will'

from nose.tools import ok_, eq_

import TreeingTools
import dendropy
from StringIO import StringIO
from itertools import ifilter

def test_generate_mrbayes_nexus():

    cmd = TreeingTools.generate_mrbayes_nexus('/path/to/alignment',
                                              '/path/to/output')
    checks = ['begin mrbayes;',
              'set autoclose=yes nowarn=yes;',
              'execute /path/to/alignment;',
              'prset aamodelpr = mixed;',
              'mcmc nchains = 3 ngen = 50000 samplefreq=1000 diagnfreq=100000 printfreq=100000 file=/path/to/output;',
              'sump;',
              'sumt;']
    for check in checks:
        yield ok_, check in cmd


def test_make_mrbayes_trees():

    seqs = [('test1', 'ATTTCTATCTATA'),
            ('test2', 'ATTTCGATCTATA'),
            ('test3', 'ATTTCGATGTATA'),
            ('test4', 'ATCTCGATGTATA'),
            ('test5', 'ATCTCGATGTATA'),
            ('test6', 'ATCTCGATGTAAA'),
            ('test7', 'ATCTCGATGTTAA'),
            ('test8', 'ATCTCGATGTTAT')]

    con_tree, all_trees = TreeingTools.make_mrbayes_trees(seqs, is_aa=False)
    for name, seq in seqs:
        node = con_tree.find_node_with_taxon_label(name)
        yield ok_, node is not None


def test_bats_format_nexus():

    taxons = dendropy.TaxonSet(['test_check%i' % i for i in range(10)])
    trop_dict = dict([('test_check%i' % i, (i % 2) == 0) for i in range(10)])
    tlist = [dendropy.treesim.uniform_pure_birth(taxons) for _ in range(20)]
    treelist = dendropy.TreeList(tlist)

    outhandle = StringIO()
    TreeingTools.bats_format_nexus(treelist, outhandle, trop_dict)

    outhandle.seek(0)
    good_lines = ifilter(lambda x: len(x.strip())>0, outhandle)
    eq_(good_lines.next().strip(), '#NEXUS', 'First line must be "#NEXUS"')
    eq_(good_lines.next().strip(), 'begin states;', 'Second line line must be "begin states;"')

    num = 0
    for num, line in enumerate(good_lines, 1):
        if line.strip() == 'End;':
            break

        tnum, state = line.strip().split()
        eq_(num, int(tnum))
        if ((num - 1) % 2) == 0:
            eq_(state, 'True')
        else:
            eq_(state, 'False')
    eq_(num, len(trop_dict) + 1, 'Some leafs were missing!')

    num = 0
    eq_(good_lines.next().strip(), 'begin trees;')
    for num, line in enumerate(good_lines, 1):
        if line.strip() == 'end;':
            break
        ok_(line.startswith('tree tree_%i' % num))
        ok_('test_check' not in line, 'Taxon names are in the tree!')
    eq_(num, len(tlist) + 1, 'Some trees were missing!')


def test_run_bats_rand_trees():

    taxons = dendropy.TaxonSet(['test_check%i' % i for i in range(100)])
    trop_dict = dict([('test_check%i' % i, (i % 2) == 0) for i in range(100)])
    tlist = [dendropy.treesim.uniform_pure_birth(taxons) for _ in range(20)]
    treelist = dendropy.TreeList(tlist)

    res = TreeingTools.run_bats(treelist, trop_dict, nreps=10)

    for row in res:
        ok_(float(row['significance']) > 0.05)

