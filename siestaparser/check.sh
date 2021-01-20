#!/bin/bash
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

python main.py test/H2O/out --annotate 2>&1 | tee test/H2O/json
python main.py test/H2O-relax/out --annotate 2>&1 | tee test/H2O-relax/json
python main.py test/Al-bulk/out --annotate 2>&1 | tee test/Al-bulk/json
python main.py test/Al-slab/out --annotate 2>&1 | tee test/Al-slab/json
python main.py test/MgO/out --annotate 2>&1 | tee test/MgO/json
python main.py test/smeagol-Au-leads/Au.out --annotate 2>&1 | tee test/smeagol-Au-leads/json
python main.py test/smeagol-Au-scregion/Au.out --annotate 2>&1 | tee test/smeagol-Au-scregion/json
