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
import sys
import logging
from glob import glob
import re

import numpy as np
from ase.data import chemical_symbols
from ase import Atoms

from nomadcore.simple_parser import (mainFunction, SimpleMatcher as SM,
                                     AncillaryParser)
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion \
    import register_userdefined_quantity, convert_unit

from siestaparser.inputvars import varlist

parser_info = {'name':'siesta-parser', 'version': '1.0'}


def siesta_energy(title, meta, **kwargs):
    return SM(r'siesta:\s*%s\s*=\s*(?P<%s__eV>\S+)' % (title, meta),
              name=meta, **kwargs)


def line_iter(fd, linepattern = re.compile(r'\s*([^#]+)')):
    # Strip off comments and whitespace, return only non-empty strings
    for line in fd:
        match = linepattern.match(line)
        if match:
            line = match.group().rstrip()
            if line:
                yield line


def get_input_metadata(inputvars_file, use_new_format):
    inputvars = {}
    blocks = {}

    varset = set(varlist)

    lower_vars = {}
    for var in varlist:
        lower_vars[var.lower()] = var

    def addvar(tokens):
        name = tokens[0]
        val = ' '.join(tokens[1:])
        name = name.lower()
        if name in lower_vars:
            name = lower_vars[name]
            inputvars[name] = val

    currentblock = None

    with open(inputvars_file) as fd:
        lines = line_iter(fd)

        for line in lines:
            tokens = line.split()
            assert len(tokens) > 0

            if tokens[0].startswith('%'):
                if tokens[0].lower() == '%block':
                    #assert currentblock == None
                    currentblock = []
                    blocks[tokens[1]] = currentblock
                elif tokens[0].lower() == '%endblock':
                    currentblock = None
                else:
                    raise ValueError('Unknown: %s' % tokens[0])
            else:
                if use_new_format:
                    if line.startswith(' '):
                        if currentblock is None:
                            continue  # Ignore.  Probably some warning
                        #assert currentblock is not None, line
                        currentblock.append(tokens)
                    else:
                        currentblock = None
                        addvar(tokens)
                else:
                    if currentblock is not None:
                        currentblock.append(tokens)
                    else:
                        addvar(tokens)

    return inputvars, blocks


"""
%block PAO.Basis                 # Define Basis set
O                     2                    # Species label, number of l-shells
 n=2   0   2                         # n, l, Nzeta
   3.305      2.479
   1.000      1.000
 n=2   1   2 P   1                   # n, l, Nzeta, Polarization, NzetaPol
   3.937      2.542
   1.000      1.000
H                     1                    # Species label, number of l-shells
 n=1   0   2 P   1                   # n, l, Nzeta, Polarization, NzetaPol
   4.709      3.760
   1.000      1.000
%endblock PAO.Basis
"""


class SiestaContext(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.fname = None  # The file that we are parsing
        self.dirname = None  # Base directory of calculations
        #self.parser = None  # The parser object
        self.format = None  # 'old' or 'new'; when parsing version

        # label and things determined by label
        self.label = None
        self.files = None  # Dict of files
        self.blocks = None  # Dict of input blocks (coords, cell, etc.)

        self.data = {}
        self.special_input_vars = {}

        self.system_meta = {}
        self.section_refs = {}  # {name: gindex, ...}
        self.simulation_type = None

    def adhoc_format_new(self, parser):
        assert self.format is None
        self.format = 'new'

    def adhoc_format_old(self, parser):
        assert self.format is None
        self.format = 'old'

    def adhoc_set_label(self, parser):
        # ASSUMPTION: the parser fIn is in the 'root' of whatever was uploaded.
        # This may not be true.  Figure out how to do this in general.
        line = parser.fIn.readline()
        assert line.startswith('reinit: System Label:')
        self.label = label = line.split()[-1]
        dirname = self.dirname

        files = {}

        for fileid in ['EIG', 'KP']:
            path = os.path.join(dirname, '%s.%s' % (label, fileid))
            if os.path.isfile(path):
                files[fileid] = path
            # Warn if files are not present?
            #
            # Also: input file
            #       input parser logfile
            #
            # what else?  We already get force/stress/positions from stdout.
        if self.format == 'new':
            inplogfiles = glob('%s/fdf-*.log' % dirname)
            assert len(inplogfiles) == 1
            if inplogfiles:
                inplogfiles.sort()
                files['inputlog'] = inplogfiles[0]
        else:
            assert self.format == 'old', self.format
            files['inputlog'] = os.path.join(dirname, 'out.fdf')
        self.files = files

    def adhoc_set_simulation_type(self, parser):
        line = parser.fIn.readline()

        if self.simulation_type is not None:
            return

        line = line.strip()

        if line.startswith('Single-point'):
            self.simulation_type = 'singlepoint'
        elif 'opt' in line or 'move' in line:
            self.simulation_type = 'optimization'
        else:
            raise ValueError('Todo: recognize simulation type "%s"' % line)

    def startedParsing(self, fname, parser):
        self.fname = fname
        path = os.path.abspath(fname)
        self.dirname = os.path.dirname(path)

    #def onClose_x_siesta_section_xc_authors(self, backend, gindex, section):

    def onClose_section_frame_sequence(self, backend, gindex, section):
        backend.addValue('frame_sequence_to_sampling_ref',
                         self.section_refs['sampling_method'])

    def onOpen_section_sampling_method(self, backend, gindex, section):
        self.section_refs['sampling_method'] = gindex

    def onOpen_section_frame_sequence(self, backend, gindex, section):
        self.section_refs['frame_sequence'] = gindex

    def onClose_section_sampling_method(self, backend, gindex, section):
        simtype = self.simulation_type
        assert simtype is not None
        if simtype == 'optimization':
            backend.addValue('sampling_method', 'geometry_optimization')
        elif simtype == 'singlepoint':
            pass
        else:
            raise ValueError('XXX: %s' % simtype)

    def onClose_section_eigenvalues(self, backend, gindex, section):
        self.read_eigenvalues(backend)

    def onOpen_section_method(self, backend, gindex, section):
        self.section_refs['method'] = gindex

    def onClose_section_method(self, backend, gindex, section):
        temp = self.special_input_vars['ElectronicTemperature']
        temp, unit = temp.split()
        assert unit == 'Ry'  # Siesta always converts to Ry here I think
        temp = float(temp)
        temp = convert_unit(temp, 'rydberg')
        backend.addValue('smearing_width', temp)

        #simtype = self.special_input_vars['MD.TypeOfRun']
        #print('SIMTYPE', simtype)
        #sdfsdf

    def onOpen_section_system(self, backend, gindex, section):
        self.section_refs['system'] = gindex

    def onClose_section_system(self, backend, gindex, section):
        data = self.data
        meta = self.system_meta

        latvec = data.pop('block_lattice_vectors', None)
        if latvec is not None:
            latvec = latvec.astype(float)
            size = self.special_input_vars['LatticeConstant']
            size, unit = size.split()
            #assert unit == 'ang', unit
            unit = {'ang': 'angstrom',
                    'Bohr': 'bohr'}[unit]

            size = float(size)
            size = convert_unit(size, unit)
            meta['simulation_cell'] = latvec * size

        cell = data.pop('outcell_ang', None)
        cell2 = data.pop('auto_unit_cell_ang', None)
        if cell2 is not None:
            cell = cell2

        if cell is not None:
            cell = cell.astype(float)
            cell = convert_unit(cell, 'angstrom')
            meta['simulation_cell'] = cell

        labels = data.pop('block_species_label', None)
        if labels is not None:
            assert labels.shape[1] == 1
            labels = labels[:, 0]
            self.labels = labels

        block_coords_and_species = self.data.pop('block_coords_and_species_from_inputlog', None)

        block_coords_and_species = data.pop('block_coords_and_species', block_coords_and_species)
        coords_and_species = data.pop('coords_and_species', None)

        if coords_and_species is None:
            coords_and_species = block_coords_and_species

        if coords_and_species is not None:
            coords = coords_and_species[:, :3].astype(float)

            unit = self.special_input_vars['AtomicCoordinatesFormat']
            if unit == 'Ang':
                coords = convert_unit(coords, 'angstrom')
            elif unit in ['Fractional', 'ScaledCartesian']:
                a = Atoms('%dX' % len(coords),
                          scaled_positions=coords,
                          cell=meta['simulation_cell'])
                coords = a.positions
            else:
                raise ValueError('Unknown: %s' % unit)
            meta['atom_positions'] = coords

            species_index = coords_and_species[:, 3].astype(int)

            atom_labels = np.array([self.labels[i - 1] for i in species_index])
            meta['atom_labels'] = atom_labels

        positions = self.data.pop('outcoord_ang', None)
        if positions is not None:
            positions = convert_unit(positions.astype(float), 'angstrom')
            meta['atom_positions'] = positions

        for key, value in meta.items():
            backend.addArrayValues(key, value)

        backend.addArrayValues('configuration_periodic_dimensions',
                               np.ones(3, bool))

        assert len(self.data) == 0, self.data

    def onClose_section_run(self, backend, gindex, section):
        pass

    def onClose_section_single_configuration_calculation(self, backend,
                                                         gindex, section):
        forces = self.data.pop('forces_ev_ang', None)
        if forces is not None:
            forces = forces.astype(float)
            forces = convert_unit(forces, 'eV/angstrom')
            backend.addArrayValues('atom_forces_free_raw', forces)

        stress = self.data.pop('stress_tensor_ev_ang', None)
        if stress is not None:
            stress = stress.astype(float)
            stress = convert_unit(stress, 'eV/angstrom^3')
            backend.addArrayValues('stress_tensor', stress)

        backend.addValue('single_configuration_to_calculation_method_ref',
                         self.section_refs['method'])
        backend.addValue('single_configuration_calculation_to_system_ref',
                         self.section_refs['system'])

    def onClose_x_siesta_section_input(self, backend, gindex, section):
        inputvars_file = self.files.get('inputlog')
        if inputvars_file is None:
            raise ValueError('no input logfile!')

        inputvars, blocks = get_input_metadata(inputvars_file,
                                               self.format == 'new')
        for varname, value in inputvars.items():
            backend.addValue('x_siesta_input_%s' % varname, value)

        for special_name in ['LatticeConstant',
                             'AtomicCoordinatesFormat',
                             'AtomicCoordinatesFormatOut',
                             'ElectronicTemperature']:
            self.special_input_vars[special_name] = inputvars.get(special_name)

        self.blocks = blocks
        self.data['block_coords_and_species_from_inputlog'] = np.array(blocks['AtomicCoordinatesAndAtomicSpecies'], str)

        authors = section['x_siesta_xc_authors']
        if authors is None:
            raise ValueError('XC authors not found!')

        assert len(authors) == 1
        authors = authors[0]

        # XXX Case sensitive?
        mapping = {'CA': ('LDA_X', 'LDA_C_PZ'),
                   'PZ': ('LDA_X', 'LDA_C_PZ'),
                   'PW92': ('LDA_X', 'LDA_C_PW'),
                   #'PW91': '',
                   'PBE': ('GGA_X_PBE', 'GGA_C_PBE'),
                   'revPBE': ('GGA_X_PBE_R', 'GGA_C_PBE'),
                   'RPBE': ('GGA_X_RPBE', 'GGA_C_PBE'),
                   #'WC': ('GGA_X_WC', ),
                   # Siesta does not mention which correlation is used with
                   # the WC functional.  Is it just the PBE one?
                   'AM05': ('GGA_X_AM05', 'GGA_C_AM05'),
                   'PBEsol': ('GGA_X_PBE_SOL', 'GGA_C_PBE_SOL'),
                   'BLYP': ('GGA_X_B88', 'GGA_C_LYP'),
                   'DF1': ('gga_x_pbe_r', 'vdw_c_df1'),
                   'DRSLL': ('gga_x_pbe_r', 'vdw_c_df1'),
                   'LMKLL': ('gga_x_rpw86', 'vdw_c_df2'),
                   'DF2': ('gga_x_rpw86', 'vdw_c_df2'),
                   'KBM': ('GGA_X_OPTB88_VDW', 'vdw_c_df1'),
                   'C09': ('GGA_X_C09X', 'vdw_c_df1'),
                   'BH': ('GGA_X_LV_RPW86', 'vdw_c_df1'),
        }
        xc = mapping.get(authors)

        if xc is None:
            raise ValueError('XC functional %s unsupported by parser'
                             % authors)

        for funcname in xc:
            gid = backend.openSection('section_XC_functionals')
            backend.addValue('XC_functional_name', funcname)
            backend.closeSection('section_XC_functionals', gid)

    def read_eigenvalues(self, backend):
        eigfile = self.files.get('EIG')
        if eigfile is None:
            return

        with open(eigfile) as fd:
            eFermi = float(next(fd).strip())
            nbands, nspins, nkpts = [int(n) for n in next(fd).split()]

            tokens = []
            for line in fd:
                tokens.extend(line.split())

        tokens = iter(tokens)

        eps = np.empty((nspins, nkpts, nbands))
        for k in range(nkpts):
            kindex = int(next(tokens))
            assert k + 1 == kindex
            for s in range(nspins):
                eps[s, k, :] = [next(tokens) for _ in range(nbands)]
        unread = list(tokens)
        assert len(unread) == 0
        assert s == nspins - 1
        assert k == nkpts - 1

        # Where does SIESTA store the occupations?
        backend.addArrayValues('eigenvalues_values', convert_unit(eps, 'eV'))

        kpfile = self.files.get('KP')
        if kpfile is None:
            return


        with open(kpfile) as fd:
            tokens = fd.read().split()

        nkpts = int(tokens[0])
        data = np.array(tokens[1:], object).reshape(-1, 5)
        coords = data[:, 1:4].astype(float)
        weights = data[:, 4].astype(float)
        backend.addArrayValues('eigenvalues_kpoints', coords)
        # XXX metadata for Fermi level?
        # k-point weights?

    def save_array(self, key, dtype=float, istart=0, iend=None,
                   unit=None):
        return get_array(key, dtype=dtype, istart=istart, iend=iend, unit=unit,
                         storage=self.data)

    def multi_sm(self, name, startpattern, linepattern, endmatcher=None,
                 conflict='fail',  # 'fail', 'keep', 'overwrite'
                 *args, **kwargs):

        pat = re.compile(linepattern)  # XXX how to get compiled pattern?
        ngroups = pat.groups

        allgroups = []
        def addline(parser):
            line = parser.fIn.readline()
            match = pat.match(line)
            assert match is not None
            thislinegroups = match.groups()
            assert len(thislinegroups) == ngroups
            allgroups.append(thislinegroups)

        def savearray(parser):
            arr = np.array(allgroups, dtype=object)
            del allgroups[:]
            if name in self.data:
                if conflict == 'fail':
                    raise ValueError('grrr %s %s' % (name, self.data[name]))
                elif conflict == 'keep':
                    return  # Do not save array
                elif conflict == 'overwrite':
                    pass
                else:
                    raise ValueError('Unknown keyword %s' % conflict)
            if arr.size > 0:
                self.data[name] = arr

        if endmatcher is None:
            endmatcher = r'.*'

        if hasattr(endmatcher, 'swapcase'):
            endmatcher = SM(endmatcher,
                            endReStr='',
                            forwardMatch=True,
                            name='%s-end' % name,
                            adHoc=savearray)

        sm = SM(startpattern,
                name=name,
                subMatchers=[
                    SM(linepattern,
                       name='%s-line' % name,
                       repeats=True,
                       forwardMatch=True,
                       required=True,
                       adHoc=addline),
                    endmatcher,
                ], **kwargs)
        return sm

context = SiestaContext()

def get_header_matcher():
    m = SM(r'',
           name='header',
           fixedStartValues={'program_name': 'Siesta',
                             'program_basis_set_type': 'numeric AOs'},
           weak=True,
           forwardMatch=True,
           subMatchers=[
               SM(r'Siesta Version: (?P<program_version>\S+)',
                  name='version',
                  adHoc=context.adhoc_format_new),
               SM(r'SIESTA\s*(?P<program_version>.+)',
                  name='alt-version', adHoc=context.adhoc_format_old),
               SM(r'Architecture\s*:\s*(?P<x_siesta_arch>.+)', name='arch'),
               SM(r'Compiler flags\s*:\s*(?P<x_siesta_compilerflags>.+)',
                  name='flags'),
           ])
    return m


def anycase(string):
    tokens = []
    for letter in list(string):
        if letter.isalpha():
            tokens.append('[%s%s]' % (letter.upper(),
                                      letter.lower()))
        else:
            tokens.append(letter)
    return ''.join(tokens)

welcome_pattern = r'\s*\*\s*WELCOME TO SIESTA\s*\*'

def get_input_matcher():
    m = SM(welcome_pattern,
           name='welcome',
           sections=['section_method', 'x_siesta_section_input'],
           fixedStartValues={'electronic_structure_method': 'DFT',
                             'smearing_kind': 'fermi'},
           subFlags=SM.SubFlags.Unordered,
           subMatchers=[
               SM(r'NumberOfAtoms\s*(?P<number_of_atoms>\d+)',
                  name='natoms'),
               context.multi_sm('block_species_label',
                                anycase(r'%block ChemicalSpeciesLabel'),
                                r'\s*\d+\s*\d+\s*(\S+)',
                            conflict='keep'),
               context.multi_sm('block_coords_and_species',
                                anycase(r'%block AtomicCoordinatesAndAtomicSpecies'),
                                r'\s*(\S+)\s*(\S+)\s*(\S+)\s*(\d+)'),
               context.multi_sm('block_lattice_vectors',
                                anycase(r'%block LatticeVectors'),
                                r'(?!%)\s*(\S+)\s*(\S+)\s*(\S+)'),
               SM(r'%s\s*(?P<x_siesta_xc_authors>\S+)' % anycase('xc\.authors'),
                  name='xc authors',
                  fixedStartValues={'x_siesta_xc_authors': 'CA'}),
               #SM(r'MD.TypeOfRun\s*(?P<x_siesta_typeofrun>\S+)',
               #   fixedStartValues={'x_siesta_typeofrun': 'none'}),
               SM(r'reinit: System Name:\s*(?P<system_name>.+)',
                  name='sysname'),
               SM(r'reinit: System Label:\s*(?P<x_siesta_system_label>\S+)',
                  name='syslabel', forwardMatch=True,
                  adHoc=context.adhoc_set_label),
               context.multi_sm('coords_and_species',
                                r'siesta: Atomic coordinates \(Bohr\) and species',
                                r'siesta:\s*(\S+)\s*(\S+)\s*(\S+)\s*(\d+)'),
               context.multi_sm('auto_unit_cell_ang',
                                r'siesta: Automatic unit cell vectors \(Ang\):',
                                r'siesta:\s*(\S+)\s*(\S+)\s*(\S+)'),
               # XXX must be an array with number of spin channels!
               #SM(r'Total number of electrons:\s*(?P<number_of_electrons>\S+)',
               #   name='nelectrons'),
       ])
    return m


step_pattern = r'\s*(Single-point calculation|Begin[^=]+=\s*\d+)'

def get_step_matcher():
    m = SM(step_pattern,
           name='step',
           forwardMatch=True,
           adHoc=context.adhoc_set_simulation_type,
           sections=['section_single_configuration_calculation'],
           subFlags=SM.SubFlags.Unordered,
           subMatchers=[
               SM(r'\s*=============+', name='====='),
               context.multi_sm('outcoord_ang',
                                r'outcoor: Atomic coordinates \(Ang\):',
                                r'\s*(\S+)\s*(\S+)\s*(\S+)'),
               context.multi_sm('outcell_ang',
                                r'outcell: Unit cell vectors \(Ang\):',
                                r'\s*(\S+)\s*(\S+)\s*(\S+)'),
               SM(r'siesta: E_KS\(eV\)\s*=\s*(?P<energy_total__eV>\S+)',
                  name='etotal'),
               context.multi_sm('forces_ev_ang',
                                r'siesta: Atomic forces \(eV/Ang\):',
                                r'(?:siesta:)?\s+\d+\s*(\S+)\s*(\S+)\s*(\S+)'),
               context.multi_sm('stress_tensor_ev_ang',
                                r'siesta: Stress tensor \(static\) \(eV/Ang\*\*3\):',
                                r'siesta:\s*(\S+)\s*(\S+)\s*(\S+)'),
               SM(r'siesta: Final energy \(eV\):',
                  endReStr='siesta:\s*Total\s*=\s*\S+',
                  name='Efinal',
                  subMatchers=[
                      siesta_energy('Band Struct\.', 'energy_sum_eigenvalues'),
                      siesta_energy('Kinetic', 'electronic_kinetic_energy'),
                      siesta_energy('Hartree', 'energy_electrostatic'),
                      #siesta_energy('Ext\. field', ''),
                      siesta_energy('Exch\.-corr\.', 'energy_XC'),
                      #siesta_energy('Ion-electron', ''),
                      #siesta_energy('Ion-Ion', ''),
                      #siesta_energy('Ekinion', ''),
                      # energy_total was matched above already
                      #siesta_energy('Total', 'energy_total'),
                      SM(r'', weak=True, name='trigger_readeig',
                         sections=['section_eigenvalues']),
                  ]),
           ])
    return m


mainFileDescription = SM(
    r'',
    name='root',
    #weak=True,
    sections=['section_run'],
    subMatchers=[
        get_header_matcher(),
        SM(r'(%s|%s)' % (welcome_pattern, step_pattern),
           name='system-section',
           #weak=True,
           forwardMatch=True,
           repeats=True,
           required=True,
           sections=['section_system', 'section_frame_sequence', 'section_sampling_method'],
           subMatchers=[
               get_input_matcher(),
               get_step_matcher(),
           ]),
        SM(r'x^',  # Make sure whole file is parsed
           name='end')
    ])


class SiestaParser():
   """ A proper class envolop for running this parser from within python. """
   def __init__(self, backend, **kwargs):
       self.backend_factory = backend

   def parse(self, mainfile):
       from unittest.mock import patch
       logging.info('siesta parser started')
       logging.getLogger('nomadcore').setLevel(logging.WARNING)
       backend = self.backend_factory("siesta.nomadmetainfo.json")
       context.reset()
       with patch.object(sys, 'argv', ['<exe>', '--uri', 'nmd://uri', mainfile]):
           mainFunction(
               mainFileDescription,
               None,
               parser_info,
               cachingLevelForMetaName={},
               superContext=context,
               superBackend=backend)

       return backend
