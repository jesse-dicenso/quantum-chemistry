#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "../libint/2e.hpp"
#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"
#include "../libdft/xc.hpp"

Matrix overlap(const std::vector<GF>& phis);
Matrix kinetic(const std::vector<GF>& phis);
Matrix nuclear(const std::vector<GF>& phis, const std::vector<int>& Zvals, const std::vector<std::vector<double>>& xyzN);
Matrix coulomb(const Matrix& P, const std::vector<std::vector<std::vector<std::vector<double>>>>& g);
Matrix fock(const Matrix& Hcore, const Matrix& J, const Matrix& K);
double nucrepl(const std::vector<int>& Z, const std::vector<std::vector<double>>& xyzN);
double E0(const XC& xc, const Matrix& Hcore, const Matrix& J);

#endif
