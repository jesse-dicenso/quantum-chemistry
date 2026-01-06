#ifndef DFT_HELPERHEADERDEF
#define DFT_HELPERHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/grid.hpp"

double density(double x, double y, double z, const Molecule& mol, const Matrix& P);
std::vector<double> density_gradient(double x, double y, double z, const Molecule& mol, const Matrix& P);

// U_VWN_c helpers
double f_zeta(double zeta);
double df_zeta(double zeta);
double VWN_alpha(double rho);
double VWN_dalpha_drho(double rho);	

#endif
