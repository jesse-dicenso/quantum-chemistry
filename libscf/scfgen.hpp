#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "../libint/2e.hpp"
#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

Matrix overlap(const std::vector<GF>& phis);
Matrix kinetic(const std::vector<GF>& phis);
Matrix nuclear(const std::vector<GF>& phis, const std::vector<int>& Zvals, const std::vector<std::vector<double>>& xyzN);
double nucrepl(const std::vector<int>& Z, const std::vector<std::vector<double>>& xyzN);

#endif
