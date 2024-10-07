#!/bin/bash

#################################
#				#
#    	  *** INPUTS ***	#
#				#
#################################

# input file name
infile="H2O.inp"

# calculation method (RHF, UHF)
method="RHF"

# basis set
basis="STO-3G"

# DIIS subspace size; if sps=0, fixed-point iterations are used
sps="0"

# convergence criterion (energy)
eps="1e-8"

# maximum number of scf iterations
max_cycles="100"

# population analysis: "lowdin", "mulliken"
pop="lowdin"

echo "Running..."
time { echo $infile; echo $method; echo $basis; echo $sps; echo $eps; echo $max_cycles; echo $pop; } | ./hartree-fock
