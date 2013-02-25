__author__ = 'will'

from nose.tools import ok_

import TreeingTools


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
