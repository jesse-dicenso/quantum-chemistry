#ifndef AUXHEADERDEF
#define AUXHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

std::vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S);
std::vector<double> Mulliken_PA(Molecule M, Matrix P, Matrix S);
void R_print_orbitals(Matrix E, Matrix C, int Nocc, int Kb);
void UR_print_orbitals(Matrix Ea, Matrix Eb, Matrix Ca, Matrix Cb, int Nocca, int Noccb, int Kb);

#endif
