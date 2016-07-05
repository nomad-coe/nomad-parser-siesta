from __future__ import print_function
import os
import sys
import setup_paths

import numpy as np
from ase.data import chemical_symbols

from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion \
    import register_userdefined_quantity, convert_unit

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

def array_matcher():

    def getarray():
        pass

    sm = SM()
    return sm

#class ArrayParser(SM):
#    def __init__(self, *args, **kwargs):
#        SM.__init__(self, *args, adHoc=get_ar**kwargs)

#def array_line_matcher(*args, **kwargs):
#    def 
#    sm = SimpleMatcher(*args, repeats=True, **kwargs)
#    return sm

def get_positions_and_labels(parser):
    rows = []
    atomic_numbers = []
    line = parser.fIn.readline()
    while line.startswith('siesta:'):
        tokens = line.split()
        rows.append([float(x) for x in tokens[1:4]])
        atomic_numbers.append(int(tokens[4]))
        line = parser.fIn.readline()
    rows = np.array(rows, float)
    labels = np.array([chemical_symbols[Z] for Z in atomic_numbers])

    b = parser.backend.superBackend
    b.addArrayValues('atom_positions', convert_unit(rows, 'bohr'))
    b.addArrayValues('atom_labels', labels)

#def get_cell(parser):
#    


def ArraySM(header, row, end, build):
    lines = []

    def addrow(parser):
        line = parser.fIn.readline()
        lines.append(line)

    def _build_array(parser):
        build(parser.backend.superBackend, lines)

    sm = SM(header,
            name='startarray',
            required=True,
            subFlags=SM.SubFlags.Sequenced,
            subMatchers=[
                SM(row, name='array', repeats=True,
                   forwardMatch=True,
                   adHoc=addrow, required=True),
                SM(end, name='endarray', required=True),
                SM(r'', adHoc=_build_array, name='dummy', forwardMatch=True)
            ])
    return sm

def build_cell(backend, lines):
    cell = np.array([[float(x) for x in line.split()] for line in lines])
    backend.addArrayValues('simulation_cell', convert_unit(cell, 'angstrom'))

def add_positions_and_labels(backend, lines):
    matrix = np.array([line.split() for line in lines], object)
    positions = matrix[:, 1:4].astype(float)
    labels = np.array([chemical_symbols[i] for i in matrix[:, 5].astype(int)])
    backend.addArrayValues('atom_positions', convert_unit(positions, 'bohr'))
    backend.addArrayValues('atom_labels', labels)

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
        ArraySM(r'siesta: Atomic coordinates \(Bohr\) and species',
                r'siesta:', r'', add_positions_and_labels),
        SM(r'\s*Single-point calculation',
           name='singleconfig',
           sections=['section_single_configuration_calculation'],
           subMatchers=[
               ArraySM(r'outcell: Unit cell vectors \(Ang\):',
                       r'\s*\S+\s*\S+\s*\S+',
                       r'\s*', build_cell),
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
