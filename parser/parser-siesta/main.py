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

from util import floating, integer

arg = sys.argv[1]
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/siesta.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath=metaInfoPath,
                                     dependencyLoader=None,
                                     extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                     uri=None)

parser_info = {'name':'siesta-parser', 'version': '1.0'}


def siesta_energy(title, meta, **kwargs):
    return SM(r'siesta:\s*%s\s*=\s*(?P<%s__eV>\S+)' % (title, meta),
              name=meta, **kwargs)


def ArraySM(header, row, build, **kwargs):

    class LineBuf:
        def __init__(self):
            self.lines = []

        def addrow(self, parser):
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
                   adHoc=linebuf.addrow, required=True),
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
    assert forces.shape[1] == 3
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


# .EIG file:
#fermi level
# nstates nspins nkpts
# kptnumber spin1eps1 spin1eps2.... spin1epsN spin2eps1 spin2eps2... spin2epsN

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
    return files

class SiestaContext(object):
    def __init__(self):
        self._is_last_configuration = False
        self.label = None
        self.fname = None
        self.dirname = None  # Base directory of calculations
        self.files = None

    def set_label(self, parser):
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

    #def is_last_configuration(self, parser):
    #    self._is_last_configuration = True

    def onClose_section_eigenvalues(self, backend, gindex, section):
        self.read_eigenvalues(backend)

    #def onClose_section_single_configuration_calculation(self, backend, gindex,
    #                                                     section):
    #    if self._is_last_configuration:
    #

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

infoFileDescription = SM(
    name='root',
    weak=True,
    startReStr='',
    # In SIESTA, the calculations are always periodic
    fixedStartValues={'program_name': 'siesta',
                      'configuration_periodic_dimensions': np.ones(3, bool)},
    sections=['section_run'],
    subFlags=SM.SubFlags.Sequenced,  # sequenced or not?
    subMatchers=[
        SM(r'Siesta Version: (?P<program_name>siesta)-(?P<program_version>\S+)',
           name='name&version', required=True),
        SM(r'reinit: System Label:\s*\S*', name='syslabel', forwardMatch=True,
           adHoc=context.set_label),
        SM(r'', weak=True, forwardMatch=True, name='system section',
           sections=['section_system'],
           subMatchers=[
               ArraySM(r'siesta: Atomic coordinates \(Bohr\) and species',
                       r'siesta:\s*\S+\s+\S+\s+\S+\s+\S+\s+\S+',
                       add_positions_and_labels),
               ArraySM(r'siesta: Automatic unit cell vectors \(Ang\):',
                       r'siesta:\s*\S+\s*\S+\s*\S+',
                       get_array('simulation_cell', float, 1, 4, unit='angstrom'))
               ]),
        SM(r'\s*Single-point calculation|\s*Begin \S+ opt\.',
           name='singleconfig',
           repeats=True,
           # XXX some of the matchers should not be in single config calculation
           sections=['section_single_configuration_calculation'],
           subFlags=SM.SubFlags.Sequenced,
           subMatchers=[
               ArraySM(r'outcoor: Atomic coordinates \(Ang\):',
                       r'\s*\S+\s*\S+\s*\S+',
                       get_array('atom_positions', float, 0, 3, unit='angstrom')),
               #ArraySM(r'outcell: Unit cell vectors \(Ang\):',
               #        r'\s*\S+\s*\S+\s*\S+',
               #        build_cell),
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
               #SM(r'outcoor: Relaxed atomic coordinates',
               #   adHoc=read_eigenvalues),
               #])
               SM(r'', weak=True, name='After singleconfigs',
                  subFlags=SM.SubFlags.Sequenced,
                  subMatchers=[
                   # There is a stupid header in the middle of nowhere which is
                   # equal to a later header, so we swallow it here:
                   #SM(r'siesta: Atomic forces \(eV/Ang\):'),
                      ArraySM(r'siesta: Atomic forces \(eV/Ang\):',
                              r'siesta:\s*[0-9]+\s+\S+\s+\S+\s+\S+',
                              get_forces, name='forces-in-single-calc'),
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


mainFunction(mainFileDescription=infoFileDescription,
             metaInfoEnv=metaInfoEnv,
             parserInfo=parser_info,
             cachingLevelForMetaName={},
             superContext=context)
