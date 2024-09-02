#!/bin/bash

#################################
#				#
#    	  *** INPUTS ***	#
#				#
#################################

# input file name and basis set
infile="CH3.inp"
basis="STO-3G"

# DIIS subspace size; if sps=0, fixed-point iterations are used
sps="0"

# convergence criterion (energy)
eps="1e-10"

# maximum number of scf iterations
max_cycles="250"

# population analysis: "lowdin", "mulliken"
pop="mulliken"

echo "Running..."
time { echo $infile; echo $basis; echo $sps; echo $eps; echo $max_cycles; echo $pop; } | ./hartree-fock
