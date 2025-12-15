#!/bin/bash

# input file name
infile=inputs/Ne.inp

# XC functional (R_, U_: HF, Slater, VWN5)
method=R_VWN5

# basis set
basis=STO-3G

# DIIS subspace size; if sps=0, fixed-point iterations are used
sps=5

# convergence criterion (energy)
eps=1e-8

# maximum number of scf iterations
max_cycles=50

# population analysis: lowdin, mulliken
pop=lowdin

echo "Running..."
time { echo $infile; echo $method; echo $basis; echo $sps; echo $eps; echo $max_cycles; echo $pop; } | ./QC-EXEC
