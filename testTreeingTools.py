__author__ = 'will'

from nose.tools import eq_, ok_

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