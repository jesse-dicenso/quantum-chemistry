#ifndef DENSITYHEADERDEF
#define DENSITYHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/quadrature.hpp"

double R_density(double x, double y, double z, const Molecule &mol, const Matrix &P);
double R_density_wrapper(double x, double y, double z, void* ctx);
double integrate_R_density(const Molecule &mol, const Matrix &P, double b_s_r, int n);

#endif
