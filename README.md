## Introduction
This is a simple quantum chemistry program that I began during Summer 2024. It is able to compute the wavefunctions of open- and closed-shell atoms and molecules in 3D, along with some other quantities of interest. It uses the McMurchie-Davidson scheme for computing all molecular integrals. This program will be continuously developed for my own learning (and perhaps for yours as well).

Most of the implementation is entirely from scratch. The only external code is from LAPACK/BLAS, specifically DSYEV (eigenvalues/eigenvectors) and DSYSV (solve linear systems, for DIIS). Everything else, including a class for matrices, was written by me!

## Current capabilities:

- Methods: RHF, UHF, DFT (RKS, UKS).

- Functionals:
  - LDA : Slater exchange, VWN5, PW92.
  - GGA : PBE exchange, PBE (restricted only).
  - MGGA: B97M-V 
  - Nonlocal: VV10

- Grid: Gauss-Chebyshev-Lebedev (100 radial points, 230 angular points).

- Basis Sets: STO-3G, def2-SVP (others may be easily added from [Basis Set Exchange](https://www.basissetexchange.org/) with some slight modifications; see libmol).

- SCF Algorithms: fixed-point, DIIS with variable subspace size.

- Population Analyses: Lowdin, Mulliken.

## How to install:
This program depends on LAPACK for some linear algebra routines (diagonalization, linear systems). If you do not already have LAPACK and BLAS installed, this can be done manually or with commands like:
```bash
sudo apt-get install libblas-dev liblapack-dev
```
DFT is also parallelized with OpenMP, which must also be installed. To build the main executable, navigate to the main/ directory, and run the following command:
```bash
make
```
The Makefile will automatically compile the program and will clean up all object files. You may need to set the compiler/flags yourself. It is important that your compiler knows where the LAPACK/BLAS/OpenMP libraries are located. The executable QC-EXEC should appear in the main/ directory, and it is ready to use.

## How to use:
The program requires a specific input file format. Sample inputs may be found in main/inputs. The input file must be located in the same directory as the executable.

The program itself is run from a bash script, run.sh. This file contains some important input settings which may be changed depending on the job. To run, simply type
```bash
bash run.sh
```
and the program will begin. It may take awhile, but it will output to the file outfile.dat. Sample output files may be found in main/sample_outputs.
