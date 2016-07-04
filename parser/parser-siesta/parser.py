from __future__ import print_function
import os
import sys
import setup_paths

from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion \
    import register_userdefined_quantity

from util import floating

arg = sys.argv[1]
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/siesta.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath=metaInfoPath,
                                     dependencyLoader=None,
                                     extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                     uri=None)

parser_info = {'name':'siesta-parser', 'version': '1.0'}

def siesta_energy(title, meta):
    return SM(r'siesta:\s*%s\s*=\s*(?P<%s__eV>\S*)' % (title, meta),
              name=meta)

infoFileDescription = SM(
    name='root',
    weak=True,
    startReStr='',
    fixedStartValues={'program_name': 'siesta'},
    sections=['section_run'],
    subFlags=SM.SubFlags.Sequenced,  # sequenced or not?
    subMatchers=[
        SM(r'Siesta Version: (?P<program_name>siesta)-(?P<program_version>\S*)',
           name='name&version'),
        SM(r'outcell: Unit cell vectors \(Ang\):', name='cell_header'),
        SM(r'\s*Single-point calculation',
           name='singleconfig',
           sections=['section_single_configuration_calculation'],
           subMatchers=[
               SM(r'siesta: Final energy \(eV\):',
                  name='energy_header',
                  subMatchers=[
                      siesta_energy('Band Struct\.', 'energy_sum_eigenvalues'),
                      siesta_energy('Kinetic', 'electronic_kinetic_energy'),
                      siesta_energy('Hartree', 'energy_electrostatic'),
                      #siesta_energy('Ext\. field', ''),
                      siesta_energy('Exch\.-corr\.', 'energy_XC'),
                      #siesta_energy('Ion-electron', ''),
                      #siesta_energy('Ion-Ion', ''),
                      #siesta_energy('Ekinion', ''),
                      siesta_energy('Total', 'energy_total')
                      ])
           ])
    ])

class SiestaContext(object):
    def startedParsing(self, *args, **kwargs):
        pass

mainFunction(mainFileDescription=infoFileDescription,
             metaInfoEnv=metaInfoEnv,
             parserInfo=parser_info,
             cachingLevelForMetaName={},
             superContext=SiestaContext())
