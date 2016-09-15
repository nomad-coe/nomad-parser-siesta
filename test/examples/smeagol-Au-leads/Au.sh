#!/bin/sh
#SBATCH -n 32
#SBATCH -t 30:00
#SBATCH -p debug
##SBATCH -U physics
##SBATCH --mem=100



mpirun -n 32 /projects/pi-ssanvito/parsons/drogheta/smeagol_acmol_soc_kshift_bond_currents/Smeagol/Src/smeagol-1.2_acmol/Src/smeagol-1.2_csg < input.fdf > Au.out
