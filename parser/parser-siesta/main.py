from __future__ import print_function
import os
import sys
import setup_paths
from glob import glob

import numpy as np
from ase.data import chemical_symbols
from ase import Atoms

from nomadcore.simple_parser import (mainFunction, SimpleMatcher as SM,
                                     AncillaryParser)
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion \
    import register_userdefined_quantity, convert_unit

from util import floating, integer
from inputvars import varlist

metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/siesta.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath=metaInfoPath,
                                     dependencyLoader=None,
                                     extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                     uri=None)

parser_info = {'name':'siesta-parser', 'version': '1.0'}


def siesta_energy(title, meta, **kwargs):
    return SM(r'siesta:\s*%s\s*=\s*(?P<%s__eV>\S+)' % (title, meta),
              name=meta, **kwargs)

import re

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

    def addvar(tokens):
        name = tokens[0]
        val = ' '.join(tokens[1:])
        if name in varset:
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
                        assert currentblock is not None
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


def ArraySM(header, row, build, **kwargs):

    class LineBuf:
        def __init__(self):
            self.lines = []

        def adhoc_addrow(self, parser):
            line = parser.fIn.readline()
            self.lines.append(line)

        def _build_array(self, parser):
            build(parser.backend.superBackend, self.lines)
            self.lines = []

    linebuf = LineBuf()
    sm = SM(header,
            name=kwargs.pop('name', 'startarray'),
            required=True,
            subFlags=SM.SubFlags.Sequenced,
            subMatchers=[
                SM(row, name='array', repeats=True,
                   forwardMatch=True,
                   adHoc=linebuf.adhoc_addrow, required=True),
                SM(r'', endReStr='', adHoc=linebuf._build_array, name='endarray',
                   forwardMatch=True)
            ],
            **kwargs)
    return sm


def errprint(func):
    def wrapper(backend, lines):
        try:
            func(backend, lines)
        except Exception:
            print('Error in %s: %d lines' % (func.__name__,
                                             len(lines)))
            for line in lines:
                print(line, file=sys.stderr)
            raise
    return wrapper

def tokenize(lines):
    return np.array([line.split() for line in lines], object)


#@errprint
#def build_cell(backend, lines):
#    cell = tokenize(lines).astype(float)
#    backend.addArrayValues('simulation_cell', convert_unit(cell, 'angstrom'))

@errprint
def get_forces(backend, lines):
    if len(lines) == 0:
        # Siesta sometimes, stupidly, writes a Forces header.
        # It then proceeds to not write any forces.
        # Thus we ignore that silliness
        return
    forces = tokenize(lines)[:, 1:].astype(float)
    assert forces.shape[1] == 3, forces.shape
    backend.addArrayValues('atom_forces', convert_unit(forces, 'eV/angstrom'))

@errprint
def get_stress(backend, lines):
    stress = tokenize(lines)[:, 1:].astype(float)
    assert stress.shape == (3, 3)
    backend.addArrayValues('stress_tensor', convert_unit(stress, 'eV/angstrom**3'))

@errprint
def get_dipole(backend, lines):
    assert len(lines) == 1
    dipole = tokenize(lines)[:, -3:].astype(float)
    # Hmm.  There is no common metadata for dipole moment for some reason.
    #print('dipole', dipole)

#def get_array(metaname, dtype, istart=0, iend=None, unit=None):
#    @errprint
#    def buildarray(backend, lines):
#        arr = tokenize(lines)
#        if iend is None:
#            arr = arr[:, istart:]
#        else:
#            arr = arr[:, istart:iend]
#        arr = arr.astype(dtype)
#        if unit is not None:
#            arr = convert_unit(arr, unit)
#        backend.addArrayValues(metaname, arr)
#    return buildarray
def get_array(metaname, dtype=float, istart=0, iend=None, unit=None,
              storage=None):
    @errprint
    def buildarray(backend, lines):
        arr = tokenize(lines)
        if iend is None:
            arr = arr[:, istart:]
        else:
            arr = arr[:, istart:iend]
        arr = arr.astype(dtype)
        if unit is not None:
            arr = convert_unit(arr, unit)
        if storage is not None:
            storage[metaname] = arr
        else:
            backend.addArrayValues(metaname, arr)
    return buildarray



@errprint
def add_positions_and_labels(backend, lines):
    matrix = tokenize(lines)
    positions = matrix[:, 1:4].astype(float)
    labels = np.array([chemical_symbols[i] for i in matrix[:, 5].astype(int)])
    backend.addArrayValues('atom_positions', convert_unit(positions, 'bohr'))
    backend.addArrayValues('atom_labels', labels)


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
        self.fname = None  # The file that we are parsing
        self.dirname = None  # Base directory of calculations
        #self.parser = None  # The parser object
        self.format = None  # 'old' or 'new'; when parsing version

        # label and things determined by label
        self.label = None
        self.files = None  # Dict of files
        self.blocks = None  # Dict of input blocks (coords, cell, etc.)

        self._is_last_configuration = False  # XXX
        self.data = {}
        self.special_input_vars = {}

        self.system_meta = {}
        self.section_refs = {}  # name: gindex

    def multi_sm0(self, pattern, var, dtypes=None, **kwargs):
        pat = re.compile(pattern)
        ngroups = pat.groups
        if dtypes is None:
            dtypes = [float] * ngroups
        else:
            assert len(dtypes) == ngroups

        class LineBuf:
            def __init__(self, ctxt):
                self.lines = []
                self.ctxt = ctxt

            def adhoc_addrow(self, parser):
                line = parser.fIn.readline()
                print('LINE', line)
                groups = pat.match(line).groups()
                results = [dtype(group)
                           for dtype, group in zip(dtypes, groups)]
                self.lines.append(results)

            def _build_array(self, parser):
                arr = np.array(self.lines)
                self.ctxt.data[var] = arr
                print('SET', var, arr)
                ksdfjksjdf
                self.lines = []

        linebuf = LineBuf(self)
        sm = SM(r'', weak=True,
                forwardMatch=True,
                subFlags=SM.SubFlags.Sequenced,
                subMatchers=[
                    SM(pattern,
                       repeats=True,
                       forwardMatch=True,
                       required=True,
                       adHoc=linebuf.adhoc_addrow),
                    SM(r'', endReStr='', adHoc=linebuf._build_array,
                       name='endarray',
                       forwardMatch=True)
                ], **kwargs)
        return sm

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
            if inplogfiles:
                inplogfiles.sort()
                files['inputlog'] = inplogfiles[-1]
        else:
            assert self.format == 'old', self.format
            files['inputlog'] = os.path.join(dirname, 'out.fdf')
        self.files = files

    def startedParsing(self, fname, parser):
        self.fname = fname
        path = os.path.abspath(fname)
        self.dirname, _ = os.path.split(path)
        #self.parser = parser

    def onClose_x_siesta_section_xc_authors(self, backend, gindex, section):
        authors = section['x_siesta_xc_authors']
        if authors is None:
            raise ValueError('XC authors not found!')

        assert len(authors) == 1
        authors = authors[0]

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
                   'BLYP': ('GGA_X_B88 + GGA_C_LYP')}
        xc = mapping.get(authors)

        if xc is None:
            raise ValueError('XC functional %s unsupported by parser'
                             % authors)

        for funcname in xc:
            gid = backend.openSection('section_XC_functionals')
            backend.addValue('XC_functional_name', funcname)
            backend.closeSection('section_XC_functionals', gid)

    def onClose_section_eigenvalues(self, backend, gindex, section):
        self.read_eigenvalues(backend)

    def onOpen_section_method(self, backend, gindex, section):
        self.section_refs['method'] = gindex

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

        block_coords_and_species = data.pop('block_coords_and_species', None)
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

        backend.addValue('single_configuration_to_calculation_method_ref',
                         self.section_refs['method'])
        backend.addValue('single_configuration_calculation_to_system_ref',
                         self.section_refs['system'])

    def onClose_x_siesta_section_input(self, backend, gindex, section):
        inputvars_file = self.files.get('inputlog')
        if inputvars_file is None:
            return

        inputvars, blocks = get_input_metadata(inputvars_file,
                                               self.format == 'new')
        for varname, value in inputvars.items():
            backend.addValue('x_siesta_input_%s' % varname, value)

        for special_name in ['LatticeConstant',
                             'AtomicCoordinatesFormat',
                             'AtomicCoordinatesFormatOut']:
            self.special_input_vars[special_name] = inputvars.get(special_name)

        self.blocks = blocks

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
               SM(r'Siesta Version: siesta-(?P<program_version>\S+)',
                  name='version',
                  adHoc=context.adhoc_format_new),
               SM(r'SIESTA\s*(?P<program_version>.+)',
                  name='alt-version', adHoc=context.adhoc_format_old),
               SM(r'Architecture\s*:\s*(?P<x_siesta_arch>.+)', name='arch'),
               SM(r'Compiler flags\s*:\s*(?P<x_siesta_compilerflags>.+)',
                  name='flags'),
           ])
    return m


welcome_pattern = r'\s*\*\s*WELCOME TO SIESTA\s*\*'

def get_input_matcher():
    m = SM(welcome_pattern,
           name='welcome',
           sections=['section_method'],
           subFlags=SM.SubFlags.Unordered,
           subMatchers=[
               SM(r'NumberOfAtoms\s*(?P<number_of_atoms>\d+)',
                  name='natoms'),
               context.multi_sm('block_species_label',
                                r'%block ChemicalSpeciesLabel',
                                r'\s*\d+\s*\d+\s*(\S+)',
                            conflict='keep'),
               context.multi_sm('block_coords_and_species',
                                r'%block AtomicCoordinatesAndAtomicSpecies',
                                r'\s*(\S+)\s*(\S+)\s*(\S+)\s*(\d+)'),
               context.multi_sm('block_lattice_vectors',
                                r'%block LatticeVectors',
                                r'(?!%)\s*(\S+)\s*(\S+)\s*(\S+)'),
               SM(r'xc.authors\s*(?P<x_siesta_xc_authors>\S+)',
                  name='xc authors',
                  fixedStartValues={'x_siesta_xc_authors': 'CA'},
                  sections=['section_method', 'x_siesta_section_xc_authors']),
               SM(r'reinit: System Label:\s*(?P<x_siesta_system_label>\S+)',
                  name='syslabel', forwardMatch=True,
                  sections=['x_siesta_section_input'],
                  adHoc=context.adhoc_set_label),
               context.multi_sm('coords_and_species',
                                r'siesta: Atomic coordinates \(Bohr\) and species',
                                r'siesta:\s*(\S+)\s*(\S+)\s*(\S+)\s*(\d+)'),
               context.multi_sm('auto_unit_cell_ang',
                                r'siesta: Automatic unit cell vectors \(Ang\):',
                                r'siesta:\s*(\S+)\s*(\S+)\s*(\S+)'),
               SM(r'Total number of electrons:\s*(?P<number_of_electrons>\S+)',
                  name='nelectrons'),
       ])
    return m


step_pattern = r'\s*(Single-point calculation|Begin[^=]+=\s*\d+)'

def get_step_matcher():
    m = SM(step_pattern,
           name='step',
           #repeats=True,
           sections=['section_single_configuration_calculation'],
           subMatchers=[
               SM(r'\s*=============+', name='====='),
               context.multi_sm('outcoord_ang',
                                r'outcoor: Atomic coordinates \(Ang\):',
                                r'\s*(\S+)\s*(\S+)\s*(\S+)'),
               context.multi_sm('outcell_ang',
                                r'outcell: Unit cell vectors \(Ang\):',
                                r'\s*(\S+)\s*(\S+)\s*(\S+)'),
               context.multi_sm('forces_ev_ang',
                                r'siesta: Atomic forces \(eV/Ang\):',
                                r'\s+\d+\s*(\S+)\s*(\S+)\s*(\S+)')
           ])
    return m


mainFileDescription = SM(
    r'',
    name='root',
    #weak=True,
    #fixedStartValues={'program_name': 'siesta',
    #                  'program_basis_set_type': 'numeric AOs',
    #                  'configuration_periodic_dimensions': np.ones(3, bool)},
    sections=['section_run'],
    subMatchers=[
        get_header_matcher(),
        SM(r'(%s|%s)' % (welcome_pattern, step_pattern),
           name='system-section',
           #weak=True,
           forwardMatch=True,
           repeats=True,
           required=True,
           sections=['section_system'],
           subMatchers=[
               get_input_matcher(),
               get_step_matcher(),
           ]),
        SM(r'x^',  # Make sure whole file is parsed
           name='end')
    ])



def main(**kwargs):
    mainFunction(mainFileDescription=mainFileDescription,
                 metaInfoEnv=metaInfoEnv,
                 parserInfo=parser_info,
                 cachingLevelForMetaName={},
                 superContext=context,
                 **kwargs)

if __name__ == '__main__':
    main()

