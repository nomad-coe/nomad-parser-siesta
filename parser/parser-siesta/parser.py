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


def siesta_energy(title, meta, **kwargs):
    return SM(r'siesta:\s*%s\s*=\s*(?P<%s__eV>\S*)' % (title, meta),
              name=meta, **kwargs)


def ArraySM(header, row, build, **kwargs):
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
                SM(r'', endReStr='', adHoc=_build_array, name='endarray',
                   forwardMatch=True)
            ],
            **kwargs)
    return sm


def tokenize(lines):
    return np.array([line.split() for line in lines], object)


def build_cell(backend, lines):
    cell = tokenize(lines).astype(float)
    backend.addArrayValues('simulation_cell', convert_unit(cell, 'angstrom'))

def get_forces(backend, lines):
    forces = tokenize(lines)[:, 2:].astype(float)
    assert forces.shape[1] == 3
    backend.addArrayValues('atom_forces', convert_unit(forces, 'eV/angstrom'))

def get_stress(backend, lines):
    stress = tokenize(lines)[:, 1:].astype(float)
    assert stress.shape == (3, 3)
    backend.addArrayValues('atom_stress', convert_unit(stress, 'eV/angstrom**3'))


def add_positions_and_labels(backend, lines):
    matrix = tokenize(lines)
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
           name='name&version', required=True),
        ArraySM(r'siesta: Atomic coordinates \(Bohr\) and species',
                r'siesta:\s*\S+\s+\S+\s+\S+\s+\S+\s+\S+',
                add_positions_and_labels),
        SM(r'\s*Single-point calculation',
           name='singleconfig',
           # XXX some of the matchers should not be in single config calculation
           sections=['section_single_configuration_calculation'],
           subFlags=SM.SubFlags.Sequenced,
           subMatchers=[
               ArraySM(r'outcell: Unit cell vectors \(Ang\):',
                       r'\s*\S+\s*\S+\s*\S+',
                       build_cell),
               SM(r'\s*scf:\s*iscf', name='scf', required=True),
               # There is a stupid header in the middle of nowhere which is
               # equal to a later header, so we swallow it here:
               SM(r'siesta: Atomic forces \(eV/Ang\):'),
               SM(r'siesta: Program\'s energy decomposition \(eV\):',
                  name='energy header 1',
                  weak=True,
                  subMatchers=[
                      siesta_energy('FreeEng', 'energy_free')
                  ]),
               SM(r'siesta: Final energy \(eV\):',
                  name='energy header 2',
                  subMatchers=[
                      siesta_energy('Band Struct\.', 'energy_sum_eigenvalues'),
                      siesta_energy('Kinetic', 'electronic_kinetic_energy'),
                      siesta_energy('Hartree', 'energy_electrostatic'),
                      #siesta_energy('Ext\. field', ''),
                      siesta_energy('Exch\.-corr\.', 'energy_XC'),
                      #siesta_energy('Ion-electron', ''),
                      #siesta_energy('Ion-Ion', ''),
                      #siesta_energy('Ekinion', ''),
                      siesta_energy('Total', 'energy_total'),
                      ]),
               ArraySM(r'siesta: Atomic forces \(eV/Ang\):',
                       r'siesta:\s*[0-9]+\s+\S+\s+\S+\s+\S+',
                       get_forces),
               ArraySM(r'siesta: Stress tensor ',#\(static\) \(eV/Ang**3\):',
                       r'siesta:\s*\S+\s+\S+\s+\S+',
                       get_stress),
               # The purpose of the following matcher is to parse all lines
               SM(r'x^', name='end')
           ])
    ])

class SiestaContext(object):
    def startedParsing(self, fname, parser):
        pass

mainFunction(mainFileDescription=infoFileDescription,
             metaInfoEnv=metaInfoEnv,
             parserInfo=parser_info,
             cachingLevelForMetaName={},
             superContext=SiestaContext())
