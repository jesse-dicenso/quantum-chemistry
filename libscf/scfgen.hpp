#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "../libint/2e.hpp"
#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

Matrix overlap(std::vector<GF> phis);
Matrix kinetic(std::vector<GF> phis);
Matrix nuclear(std::vector<GF> phis, std::vector<int> Zvals, std::vector<std::vector<double>> xyzN);
double nucrepl(std::vector<int> Z, std::vector<std::vector<double>> xyzN);

#endif
