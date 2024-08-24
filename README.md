## Introduction
This is a simple Hartree-Fock program that I wrote during Summer 2024. It is able to compute the RHF wavefunction of closed-shell atoms and molecules in 3D, along with some other quantities of interest. It uses the McMurchie-Davidson scheme for computing all molecular integrals.

Most of the implementation is entirely from scratch. The only external code is from LAPACK/BLAS, specifically for DSYEV (eigenvalues/eigenvectors) and DSYSV (solve linear systems). Everything else, including a class for matrices, was written by me!

## Current capabilities:

-Basis Sets: STO-3G

-SCF Algorithms: fixed-point, DIIS

-Population Analyses: Lowdin, Mulliken

Currently, the basis sets are only implemented for H, He, C, and O; I am still working on implementing a system so basis sets may be read in from files.

## How to install:
This program depends on LAPACK for some linear algebra routines (diagonalization, linear systems). If you do not already have LAPACK and BLAS installed, this can be done manually or with commands like:
```bash
sudo apt-get install libblas-dev liblapack-dev
```

With LAPACK and BLAS installed, simply download all of the program directories/files, navigate to the /hartree-fock/ directory, and run the following command:
```bash
make hartree-fock
```

The Makefile will automatically compile the program and will clean up all object files. The program is compiled with g++; it is important that your compiler knows where the LAPACK/BLAS libraries are located. The executable hartree-fock should appear in the /hartree-fock/ directory, and it is ready to use.

## How to use:
The program requires a specific input file format. Sample inputs may be found in /hartree-fock/sample_inputs. The input file must be located in the same directory as the executable.

The program itself is run from a bash script, run.sh. This file contains some important input settings which may be changed depending on the job. To run, simply type
```bash
bash run.sh
```
and the program will begin. It may take awhile, but it will output to the file outfile.dat. Sample output files may be found in the directory /hartree-fock/sample_outputs.

There may be some numerical instabilities, but it generally works well for small molecules.
