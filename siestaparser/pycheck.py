# Copyright 2016-2018 Ask Hjorth Larsen, Fawzi Mohamed
# 
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#
# Main author and maintainer: Ask Hjorth Larsen <asklarsen@gmail.com>

from __future__ import print_function
import os

from main import main

testdir = '../../test/examples'
dirnames = ['H2O',
            'H2O-relax',
            'Al-slab',
            'Al-uc',
            'Fe',
            'MgO']


for dirname in dirnames:
    fname = os.path.join(testdir, dirname, 'out')
    with open('out.pycheck.%s.txt' % dirname, 'w') as outfd:
        main(mainFile=fname, outF=outfd)
