#ifndef DENSITYHEADERDEF
#define DENSITYHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/quadrature.hpp"

// Functions to integrate
double R_atomic_density(double x, double y, double z, const Molecule &mol, const Matrix &P);

// Wrappers for above functions (use these as inputs to integrate_R)
double R_atomic_density_wrapper(double x, double y, double z, void* ctx);

// Integration Routines for above functions
double integrate_R(const Molecule &mol, const Matrix &P, int n);
double becke_weight(double x, double y, double z, int l, int k, const Molecule &mol);
double becke_step(double mu, int k);

extern const double bragg_slater_radii[118];

#endif
