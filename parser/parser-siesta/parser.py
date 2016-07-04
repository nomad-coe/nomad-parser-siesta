from __future__ import print_function
import sys
import setup_paths

from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion \
    import register_userdefined_quantity


arg = sys.argv[1]

infoFileDescription = SM(
    name='root',
    weak=True,
    startReStr='',
    fixedStartValues={'program_name': 'siesta'},
    sections=['section_run'],
    subFlags=SM.SubFlags.Sequenced,
    subMatchers=[
        SM(r'Siesta Version: (?P<program_version>\S*)')
    ])


