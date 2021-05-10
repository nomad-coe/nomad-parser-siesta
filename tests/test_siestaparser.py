#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest

from nomad.datamodel import EntryArchive
from siestaparser import SiestaParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return SiestaParser()


def test_single_point(parser):
    archive = EntryArchive()

    parser.parse('tests/data/Fe/out', archive, None)

    sec_run = archive.section_run[0]
    assert sec_run.program_version == 'siesta-4.0--500'

    sec_system = archive.section_run[0].section_system[0]
    assert sec_system.lattice_vectors[2][2].magnitude == approx(-1.435e-10)
    assert sec_system.atom_labels[0] == 'Fe'

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert sec_scc.energy_total.magnitude == approx(-1.25299539e-16)


def test_realax(parser):
    archive = EntryArchive()

    parser.parse('tests/data/H2O-relax/out', archive, None)

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 7
    assert sec_systems[3].lattice_vectors[1][1].magnitude == approx(5.692925e-10)
    assert sec_systems[6].atom_positions[2][0].magnitude == approx(-7.6940893e-11)
    assert sec_systems[0].atom_labels[0] == 'O'

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 7
    assert sec_sccs[5].energy_method_current.magnitude == approx(-7.48200466e-17)
    assert sec_sccs[1].energy_total.magnitude == approx(-7.48202036e-17)
    assert sec_sccs[2].atom_forces[1][1].magnitude == approx(-2.06092787e-10)
