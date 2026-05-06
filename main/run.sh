#!/bin/bash

# input/output file names
infile=inputs/acetaldehyde.inp
outfile=outfile.dat

# XC functional (R_, U_: HF, HFSNX, Slater, VWN5, PW92, PBE_X, PBE (R_ only), B97M-V (U_ only))
method=R_PBE0-DH

# basis set
basis=STO-3G

# DIIS subspace size; if sps=0, fixed-point iterations are used
sps=5

# convergence criterion (energy)
eps=1e-8

# maximum number of scf iterations
max_cycles=50

# population analysis: mulliken, lowdin
pop=lowdin

# screening threshold for ERIs and SNX (if applicable)
eri_thresh=1e-10
snx_thresh_e=1e-10
snx_thresh_k=1e-8

export OMP_NUM_THREADS=$(nproc);
#export OMP_NUM_THREADS=1;
{ time ./QC-EXEC $infile $method $basis $sps $eps $max_cycles $pop $eri_thresh $snx_thresh_e $snx_thresh_k; } 2>&1 | tee $outfile
