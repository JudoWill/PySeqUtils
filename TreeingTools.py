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


@contextlib.contextmanager
def tmp_directory(*args, **kwargs):
    """A context manager which changes the working directory to the given
    path, and then changes it back to its previous value on exit.

    """
    path = mkdtemp(*args, **kwargs)
    try:
        yield path + '/'
    finally:
        shutil.rmtree(path)


def make_mrbayes_trees(input_seqs, mrbayes_kwargs=None, is_aa=True):
    """Takes an ALIGNED set of sequences and generates a phylogenetic tree using MrBayes.
    """

    with tmp_directory() as tmpdir:
        align_file = tmpdir + 'seqalign.nxs'
        mrbayes_cmd_file = tmpdir + 'analysis.nxs'
        multi_prob = tmpdir + 'seqalign.nxs.trprobs'
        cons_file = tmpdir + 'seqalign.nxs.con.tre'

        with open(align_file, 'w') as handle:
            write_nexus_alignment(input_seqs, handle, is_aa=is_aa)

        with open(mrbayes_cmd_file, 'w') as handle:
            if mrbayes_kwargs is None:
                mrbayes_kwargs = {}
            txt = generate_mrbayes_nexus(align_file, align_file, is_aa=is_aa, **mrbayes_kwargs)
            handle.write(txt)

        cmd = '/home/will/mb ' + mrbayes_cmd_file
        check_output(shlex.split(cmd))

        with open(multi_prob) as handle:
            trees = dendropy.TreeList.get_from_stream(handle, schema='nexus')
        with open(cons_file) as handle:
            con_tree = dendropy.Tree.get_from_stream(handle, schema='nexus')

    return con_tree, trees


def generate_mrbayes_nexus(alignment_path, output_path,
                           nchains=3, ngen=50000, samplefreq=1000,
                           is_aa=True):
    """Generates the NEXUS command to control MrBayes in the form that I usually use. This will likely be expanded as I
     have to include more issues.
    """

    cmd = """begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute %(align)s;
   %(model)s;
   mcmc nchains = %(nchains)i ngen = %(ngen)i samplefreq=%(samplefreq)i diagnfreq=100000 printfreq=100000 file=%(out)s;
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


def run_bats(treeset, trop_dict, nreps=5000):
    """Runs the BatS analysis on a treeset.
    """

    nstates = len(set(trop_dict.values()))
    with NTF() as handle:
        bats_format_nexus(treeset, handle, trop_dict)
        handle.flush()
        os.fsync(handle)
        cmd = 'java -jar /home/will/BaTS_beta_build2.jar single %s %i %i'
        out = check_output(shlex.split(cmd % (handle.name, nreps, nstates)))

    handle = StringIO(out)
    headers = []
    for line in handle:
        if line.startswith('Stat'):
            headers = line.strip().split('\t')
            break
    return list(csv.DictReader(handle, fieldnames=headers, delimiter='\t'))[:-2]
