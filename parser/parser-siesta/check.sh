#!/bin/bash

python main.py test/H2O/out --annotate
python main.py test/H2O-relax/out --annotate
python main.py test/Al-bulk/out --annotate
python main.py test/Al-slab/out --annotate
python main.py test/MgO/out --annotate
python main.py test/smeagol-Au-leads/Au.out --annotate
python main.py test/smeagol-Au-scregion/Au.out --annotate
