__author__ = 'will'

import dendropy


def make_tree(input_seqs):
    """Takes an ALIGNED set of sequences and generates a phylogenetic tree using MrBayes.
    """

    pass


def generate_mrbayes_nexus(alignment_path, output_path,
                           nchains=3, ngen=50000, samplefreq=1000,
                           is_aa=True):

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