#ifndef DFT_HELPERHEADERDEF
#define DFT_HELPERHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/grid.hpp"

struct density_context{
	const Molecule 	*molecule;
	const Matrix	*Pmatrix;
	int 			idx1;
	int 			idx2;
};

// Functions to integrate
double density(double x, double y, double z, void* ctx);
std::vector<double> density_gradient(double x, double y, double z, void* ctx);

double R_Slater_X_integrand(double x, double y, double z, void* ctx);
double U_Slater_X_integrand(double x, double y, double z, void* ctx);

// Integration Routines (mostly for testing)
double integrate_density(const grid &g, const Molecule &mol, const Matrix &P);

#endif
