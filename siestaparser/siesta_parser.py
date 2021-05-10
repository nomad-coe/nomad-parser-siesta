#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from nomad.units import ureg
from nomad.parsing.file_parser import BasicParser


class SiestaParser(BasicParser):
    def __init__(self):
        re_f = r'\-*\d+\.\d+E*e*\-*\+*\d*'

        super().__init__(
            specifications=dict(
                name='parsers/siesta', code_name='Siesta', code_homepage='https://departments.icmab.es/leem/siesta/',
                mainfile_contents_re=(
                    r'(Siesta Version: siesta-|SIESTA [0-9]\.[0-9]\.[0-9])|'
                    r'(\*\s*WELCOME TO SIESTA\s*\*)')),
            units_mapping=dict(energy=ureg.eV, length=ureg.angstrom),
            program_version=r'Siesta Version\: (siesta\S+)',
            lattice_vectors=r'outcell\: Unit cell vectors \(Ang\)\:([\s\S]+?)\n *\n',
            # TODO support for various input, output formats
            atom_labels_atom_positions_scaled=r'outcoor\: Atomic coordinates \((?:fractional|scaled)\)\:([\s\S]+?\n *\n)',
            atom_labels_atom_positions=r'outcoor\: Atomic coordinates \(Ang\)\:([\s\S]+?\n *\n)',
            energy_method_current=rf'siesta\: E_KS\(eV\) = +({re_f})',
            energy_total=rf'siesta\: Etot +\= +({re_f})',
            atom_forces=r'siesta: Atomic forces \(eV/Ang\)\:([\s\S]+?)\-{10}')
