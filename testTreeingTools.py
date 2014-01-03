__author__ = 'will'

from nose.tools import ok_, eq_

from Bio.Alphabet import generic_dna, generic_protein
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
              'sump;',
              'sumt;']
    for check in checks:
        yield ok_, check in cmd, 'Missing: "%s"' % check


def tree_seqs():

    seqs = [('test1',   'ATTTCTATCTATA'),
            ('test1-1', 'ATTTtTATtTATA'),
            ('test1-2', 'ATcTtTATtTATA'),
            ('test1-3', 'ATTctTATtTATA'),
            ('test1-4', 'ATTcCTATCTATA'),
            ('test2',   'TTTTCCCGGTGTG'),
            ('test2-1', 'TTTTCCCGtTGgG'),
            ('test2-2', 'TTTTCgCGtTGgG'),
            ('test2-3', 'TTTTCCgGtTGgG'),
            ('test2-4', 'TTTTCCgGtTcTG'),
            ('test3',   'AATGATCGATTTA'),
            ('test3-1', 'AcTGATgGATTTA'),
            ('test3-2', 'AcTGAggGATTTA'),
            ('test3-3', 'AcTGATgGATgTA'),
            ('test3-4', 'AATGtTgGAgTTA')]
    return seqs


def check_tree(tree):

    with open('/tmp/ascii_tree.txt', 'w') as handle:
        handle.write(tree.as_ascii_plot())

    for name in ['test1', 'test2', 'test3']:
        node = tree.find_node_with_taxon_label(name)
        children = [node for node in node.sister_nodes() if node.taxon]
        children += [node for node in node.child_nodes() if node.taxon]
        #yield ok_, len(children) > 0, '%s did not have the right number of children: %i' % (name, len(children))
        for child in children:
            if child.taxon:
                yield ok_, child.taxon.label.startswith(name), 'The tree is wrong!'


def test_make_mrbayes_trees():

    seqs = tree_seqs()
    con_tree, all_trees = TreeingTools.make_mrbayes_trees(seqs, is_aa=False)
    for tst in check_tree(con_tree):
        #yield tst
        pass


def test_phylip_tree():

    seqs = tree_seqs()
    tree, _ = TreeingTools.phylip_tree(seqs, alphabet=generic_dna)
    for tst in check_tree(tree):
        yield tst


def test_fast_tree():

    seqs = tree_seqs()
    tree = TreeingTools.run_FastTree(seqs, alphabet=generic_dna)
    for tst in check_tree(tree):
        yield tst


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


