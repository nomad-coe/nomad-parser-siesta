from __future__ import print_function
import os

from main import main

testdir = '../../test/examples'
dirnames = ['H2O',
            'h2o_relax',
            'Al_slab',
            'Al_uc',
            'Fe',
            'MgO']


for dirname in dirnames:
    fname = os.path.join(testdir, dirname, 'out')
    with open('out.pycheck.%s.txt' % dirname, 'w') as outfd:
        main(mainFile=fname, outF=outfd)
