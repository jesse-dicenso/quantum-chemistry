#!/bin/bash

# input file name
infile=inputs/H2O.inp

# XC functional (R_, U_: HF, HFSNX, Slater, VWN5, PW92, PBE_X, PBE (R_ only), B97M-V (U_ only))
method=R_HFSNX

# basis set
basis=STO-3G

# DIIS subspace size; if sps=0, fixed-point iterations are used
sps=5

# convergence criterion (energy)
eps=1e-8

# screening threshold for ERIs
int_thresh=1e-12

# maximum number of scf iterations
max_cycles=50

# population analysis: mulliken, lowdin
pop=lowdin

export OMP_NUM_THREADS=$(nproc);
{ time ./QC-EXEC $infile $method $basis $sps $eps $int_thresh $max_cycles $pop; } 2>&1 | tee outfile.dat
