#!/bin/bash

python main.py test/H2O/out --annotate > test/H2O/json 2>&1
python main.py test/H2O-relax/out --annotate > test/H2O-relax/json 2>&1
python main.py test/Al-bulk/out --annotate > test/Al-bulk/json 2>&1
python main.py test/Al-slab/out --annotate > test/Al-slab/json 2>&1
python main.py test/MgO/out --annotate > test/MgO/json 2>&1
python main.py test/smeagol-Au-leads/Au.out --annotate > test/smeagol-Au-leads/json 2>&1
python main.py test/smeagol-Au-scregion/Au.out --annotate > test/smeagol-Au-scregion/json 2>&1
