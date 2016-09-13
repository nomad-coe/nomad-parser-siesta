from __future__ import print_function
import os
import sys
import setup_paths
from glob import glob

import numpy as np
from ase.data import chemical_symbols

from nomadcore.simple_parser import (mainFunction, SimpleMatcher as SM,
                                     AncillaryParser)
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion \
    import register_userdefined_quantity, convert_unit

from util import floating, integer

metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/siesta.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath=metaInfoPath,
                                     dependencyLoader=None,
                                     extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                     uri=None)

parser_info = {'name':'siesta-parser', 'version': '1.0'}


def siesta_energy(title, meta, **kwargs):
    return SM(r'siesta:\s*%s\s*=\s*(?P<%s__eV>\S+)' % (title, meta),
              name=meta, **kwargs)


def get_input_metadata(inputvars_file):
    inputvars = {}
    with open(inputvars_file) as fd:
        for line in fd:
            if not line or line.startswith('%') or line[0].isspace():
                continue

            line = line.split('#', 1)[0]  # Strip off comment
            tokens = line.split()
            if not tokens:
                continue
            tokens = line.split()
            if not tokens:
                continue

            metaname = 'x_siesta_input_%s' % tokens[0]
            value = ' '.join(tokens[1:])
            inputvars[metaname] = value
    return inputvars


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

def get_array(metaname, dtype, istart=0, iend=None, unit=None):
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


def locate_files(dirname, label):
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
    inplogfiles = glob('%s/fdf-*.log' % dirname)
    if inplogfiles:
        inplogfiles.sort()
        files['inputlog'] = inplogfiles[-1]
    return files

class SiestaContext(object):
    def __init__(self):
        self._is_last_configuration = False
        self.label = None
        self.fname = None
        self.dirname = None  # Base directory of calculations
        self.files = None
        self.parser = None

    def adhoc_set_label(self, parser):
        # ASSUMPTION: the parser fIn is in the 'root' of whatever was uploaded.
        # This may not be true.  Figure out how to do this in general.
        line = parser.fIn.readline()
        assert line.startswith('reinit: System Label:')
        self.label = line.split()[-1]
        self.files = locate_files(self.dirname, self.label)

    def startedParsing(self, fname, parser):
        self.fname = fname
        path = os.path.abspath(fname)
        self.dirname, _ = os.path.split(path)
        self.parser = parser

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
                   'PBE': ('LDA_X_PBE', 'LDA_C_PBE'),
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

    def onClose_x_siesta_section_input(self, backend, gindex, section):
        inputvars_file = self.files.get('inputlog')
        if inputvars_file is None:
            return

        inputvars = get_input_metadata(inputvars_file)
        for metaname, value in inputvars.items():
            backend.addValue(metaname, value)

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


context = SiestaContext()

from inputvars import inputvars

input_vars_description = SM(
    name='root',
    weak=True,
    startReStr='',
    subMatchers= [
        SM(r'%{0}\s+(?P<x_siesta_input_{0}>\S+)'.format(inputvar),
           name=inputvar)
        for inputvar in inputvars]
    )

infoFileDescription = SM(
    name='root',
    weak=True,
    startReStr='',
    # In SIESTA, the calculations are always periodic
    fixedStartValues={'program_name': 'siesta',
                      'configuration_periodic_dimensions': np.ones(3, bool)},
    sections=['section_run', 'x_siesta_section_input'],
    subFlags=SM.SubFlags.Sequenced,  # sequenced or not?
    subMatchers=[
        SM(r'Siesta Version: siesta-(?P<program_version>\S+)',
           name='name&version', required=True),
        SM(r'xc.authors\s*(?P<x_siesta_xc_authors>\S+)',
           name='xc authors',
           fixedStartValues={'x_siesta_xc_authors': 'CA'},
           sections=['section_method', 'x_siesta_section_xc_authors']),
        SM(r'reinit: System Label:\s*(?P<x_siesta_system_label>\S+)', name='syslabel', forwardMatch=True,
           adHoc=context.adhoc_set_label),
        SM(r'\s*Single-point calculation|\s*Begin \S+ opt\.',
           name='singleconfig',
           repeats=True,
           # XXX some of the matchers should not be in single config calculation
           sections=['section_single_configuration_calculation', 'section_system'],
           subFlags=SM.SubFlags.Sequenced,
           subMatchers=[
               SM(r'', weak=True, forwardMatch=True, name='system section',
                  subMatchers=[
                      ArraySM(r'siesta: Atomic coordinates \(Bohr\) and species',
                              r'siesta:\s*\S+\s+\S+\s+\S+\s+\S+\s+\S+',
                              add_positions_and_labels),
                      ArraySM(r'siesta: Automatic unit cell vectors \(Ang\):',
                              r'siesta:\s*\S+\s*\S+\s*\S+',
                              get_array('simulation_cell', float, 1, 4, unit='angstrom'))
                  ]),
               ArraySM(r'outcoor: Atomic coordinates \(Ang\):',
                       r'\s*\S+\s*\S+\s*\S+',
                       get_array('atom_positions', float, 0, 3, unit='angstrom')),
               SM(r'\s*scf:\s*iscf', name='scf', required=True),
               SM(r'SCF cycle converged after\s*%s\s*iterations'
                  % integer('number_of_scf_iterations'), name='scf-iterations'),
               SM(r'siesta: E_KS\(eV\)\s*=\s*%s' % floating('energy_free__eV'),
                  name='E_KS'),
               ArraySM(r'siesta: Atomic forces \(eV/Ang\):',
                       r'\s*[0-9]+\s+\S+\s+\S+\s+\S+',
                       get_forces, name='forces-in-opt-step'),
               SM(r'Target enthalpy', name='end-singleconfig'),
               # As we get past the last singleconfig, we need to add all the stuff
               # that pertains only to the last one, basically eigs and occs.
               # Therefore we read a bit past the last config and try to trigger
               # this as appropriate.
               SM(r'', weak=True, name='After singleconfigs',
                  subFlags=SM.SubFlags.Sequenced,
                  subMatchers=[
                   # There is a stupid header in the middle of nowhere which is
                   # equal to a later header, so we swallow it here:
                   #SM(r'siesta: Atomic forces \(eV/Ang\):'),
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
                             SM(r'', weak=True, name='trigger_readeig',
                                sections=['section_eigenvalues']),
                         ]),
                      ArraySM(r'siesta: Atomic forces \(eV/Ang\):',
                              r'siesta:\s*[0-9]+\s+\S+\s+\S+\s+\S+',
                              get_array('atom_forces', float, 2, 5,
                                        unit='eV/angstrom'),
                              name='forces-in-single-calc'),
                      ArraySM(r'siesta: Stress tensor \(static\) \(eV/Ang\*\*3\):',
                              r'siesta:\s*\S+\s+\S+\s+\S+',
                              get_stress),
                      ArraySM(r'siesta: Electric dipole \(Debye\)\s*=',
                              r'siesta: Electric dipole \(Debye\)\s*=',
                              get_dipole,
                              forwardMatch=True),
                      # The purpose of the following matcher is to parse all lines
                      SM(r'x^', name='end')
                  ])
           ])
    ])

def main(**kwargs):
    mainFunction(mainFileDescription=infoFileDescription,
                 metaInfoEnv=metaInfoEnv,
                 parserInfo=parser_info,
                 cachingLevelForMetaName={},
                 superContext=context,
                 **kwargs)

if __name__ == '__main__':
    main()

