#!/bin/bash

# Input file
infile="CH4.inp"
# SCF Loop
eps="1e-6"
max_cycles="50"
# Population Analysis: "lowdin", "mulliken"
pop="lowdin"

time { echo $infile; echo $eps; echo $max_cycles; echo $pop; } | ./hartree-fock
