#ifndef DENSITYHEADERDEF
#define DENSITYHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/grid.hpp"

// Functions to integrate
double R_density(double x, double y, double z, void* ctx);

// Integration Routines for above functions
double integrate_R_density(const grid &g, const Molecule &mol, const Matrix &P);

#endif
