#ifndef DENSITYHEADERDEF
#define DENSITYHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/grid.hpp"

// Functions to integrate
double density(double x, double y, double z, void* ctx);
std::vector<double> density_gradient(double x, double y, double z, void* ctx);

// Integration Routines for above functions
double integrate_density(const grid &g, const Molecule &mol, const Matrix &P);

#endif
