#ifndef DFT_HELPERHEADERDEF
#define DFT_HELPERHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/grid.hpp"

// Useful helper functions
double density(double x, double y, double z, const Molecule& mol, const Matrix& P);
std::vector<double> density_gradient(double x, double y, double z, const Molecule& mol, const Matrix& P);

#endif
