#!/bin/bash

#####################################
#				    #
#    MAIN INPUTS TO hartree-fock    #
#				    #
#####################################

# input file name
infile="HeH+.inp"

# DIIS subspace size; if sps=0, fixed-point iterations are used
sps=0

# convergence criterion (energy)
eps="1e-8"

# maximum number of scf iterations
max_cycles="250"

# population analysis: "lowdin", "mulliken"
pop="lowdin"

echo "Please be patient; this may take awhile!"
time { echo $infile; echo $sps; echo $eps; echo $max_cycles; echo $pop; } | ./hartree-fock
