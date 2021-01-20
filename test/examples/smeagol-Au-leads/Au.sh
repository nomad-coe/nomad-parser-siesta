#!/bin/sh
##
## Copyright The NOMAD Authors.
##
## This file is part of NOMAD.
## See https://nomad-lab.eu for further info.
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##

#SBATCH -n 32
#SBATCH -t 30:00
#SBATCH -p debug
##SBATCH -U physics
##SBATCH --mem=100



mpirun -n 32 /projects/pi-ssanvito/parsons/drogheta/smeagol_acmol_soc_kshift_bond_currents/Smeagol/Src/smeagol-1.2_acmol/Src/smeagol-1.2_csg < input.fdf > Au.out
