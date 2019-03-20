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

# TEMPLATES TAKEN FROM OCTOPUS.  KEEP SYNCHRONIZED.
json_header = """{
  "type": "nomad_meta_info_1_0",
  "description": "autogenerated nomad meta info for %(parser)s parser.  The file which generates this maintained in the parser's git repository.",
  "dependencies": [ {
      "relativePath": "common.nomadmetainfo.json"
    }],
  "metaInfos": [ %(info)s]
}"""


json_template = """{
      "description": "%(description)s",
      "dtypeStr": "%(dtypeStr)s",
      "name": "%(name)s",
      "repeats": false,
      "shape": [],
      "superNames": [
        "%(supername)s"
      ]
    }"""

json_section_template = """{
      "description": "%(description)s",
      "kindStr": "type_abstract_document_content",
      "name": "%(name)s",
      "superNames": [
        "x_siesta_section_input"
      ]
    }, {
      "description": "input section",
      "kindStr": "type_section",
      "name": "x_siesta_section_input",
      "superNames": [
        "section_method"
      ]
    }"""
# END OF TEMPLATES FROM OCTOPUS

import sys
#fname = sys.argv[1]

inputvars_fd = open('inputvars.py', 'w')
print('varlist = [', file=inputvars_fd)

def getlines():
    fd = open('test/H2O/fdf-95611.log')
    for line in fd:
        if not line[:1].isalpha():  # Get rid of blocks
            continue
        yield line

    withinblock = False
    fd = open('test/smeagol-Au-leads/out.fdf')
    for line in fd:
        if line.startswith('%block'):
            withinblock = True
        elif line.startswith('%endblock'):
            withinblock = False
        if withinblock:
            continue
        yield line

seen = set()

#with open(fname) as fd:
#for line in getlines():
if 1:
    varnames = []
    for line in getlines():

        # Get rid of comments
        line = line.strip().rsplit('#', 1)[0].strip()
        if not line:
            continue
        varname = line.split()[0]
        lowervarname = varname.lower()
        if lowervarname in seen:
            continue  # Already registered

        seen.add(lowervarname)
        varnames.append(varname)
    varnames = list(set(varnames))
    varnames.sort()

    jsontokens = []
    for var in varnames:
        print("    '%s'," % var, file=inputvars_fd)
        json = json_template % dict(description=r'siesta input variable \"%s\"'
                                    % var,
                                    name='x_siesta_input_%s' % var,
                                    dtypeStr='C',
                                    supername='x_siesta_input')
        jsontokens.append(json)
    section = json_section_template % dict(description='siesta input '
                                           'variables',
                                           name='x_siesta_input')
    jsontokens.append(section)

    txt = json_header % dict(parser='siesta',
                             info=', '.join(jsontokens))
    print(txt)

print(']', file=inputvars_fd)