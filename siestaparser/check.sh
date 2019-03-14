#!/bin/bash

python main.py test/H2O/out --annotate 2>&1 | tee test/H2O/json
python main.py test/H2O-relax/out --annotate 2>&1 | tee test/H2O-relax/json
python main.py test/Al-bulk/out --annotate 2>&1 | tee test/Al-bulk/json
python main.py test/Al-slab/out --annotate 2>&1 | tee test/Al-slab/json
python main.py test/MgO/out --annotate 2>&1 | tee test/MgO/json
python main.py test/smeagol-Au-leads/Au.out --annotate 2>&1 | tee test/smeagol-Au-leads/json
python main.py test/smeagol-Au-scregion/Au.out --annotate 2>&1 | tee test/smeagol-Au-scregion/json
