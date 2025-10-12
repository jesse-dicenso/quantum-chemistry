#!/bin/bash

# input file name
infile="H2O.inp"

# calculation method (RHF, UHF)
method="RHF"

# basis set
basis="def2-SVP"

# DIIS subspace size; if sps=0, fixed-point iterations are used
sps="3"

# convergence criterion (energy)
eps="1e-6"

# maximum number of scf iterations
max_cycles="100"

# population analysis: "lowdin", "mulliken"
pop="lowdin"

echo "Running..."
time { echo $infile; echo $method; echo $basis; echo $sps; echo $eps; echo $max_cycles; echo $pop; } | ./hartree-fock
