#ifndef AUXHEADERDEF
#define AUXHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

std::vector<double> Lowdin_PA  (const Molecule& M, const Matrix& P, const Matrix& S);
std::vector<double> Mulliken_PA(const Molecule& M, const Matrix& P, const Matrix& S);

void R_print_orbital_energies(const Matrix& E, int Nocc, int Kb);
void R_print_orbitals(const Matrix& C, int Nocc, int Kb);
void UR_print_orbitals(const Matrix& Ea, const Matrix& Eb, int Nocca, int Noccb, int Kb);
void UR_print_orbital_energies(const Matrix& Ca, const Matrix& Cb, int Nocca, int Noccb, int Kb);

#endif
