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
        "section_run"
      ]
    }"""
# END OF TEMPLATES FROM OCTOPUS

import sys
fname = sys.argv[1]

with open(fname) as fd:
    varnames = []
    for line in fd:
        if not line[:1].isalpha():  # Get rid of blocks
            continue

        # Get rid of comments
        line = line.strip().rsplit('#', 1)[0]
        varname = line.split()[0]
        varnames.append(varname)
    varnames = list(set(varnames))
    varnames.sort()

    jsontokens = []
    for var in varnames:
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
