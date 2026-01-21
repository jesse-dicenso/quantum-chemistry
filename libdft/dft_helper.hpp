#ifndef DFT_HELPERHEADERDEF
#define DFT_HELPERHEADERDEF

#include "../libmath/linalg.hpp"
#include "../libgf/gf.hpp"
#include "../libmol/mol.hpp"
#include "../libgrid/grid.hpp"

double density(const std::vector<double>& phis, const Matrix& P);
std::vector<double> density_gradient(const std::vector<double>& phis, const std::vector<double>& gpx, const std::vector<double>& gpy, 
									 const std::vector<double>& gpz, const Matrix& P);
double ke_density(const std::vector<double>& gpx, const std::vector<double>& gpy, const std::vector<double>& gpz, const Matrix& P);

// U_LDA_c helpers
double f_zeta(double zeta);
double df_zeta(double zeta);
double VWN_alpha(double x);
double VWN_dalpha_drho(double x, double n);	
double PW92_alpha(double rs);
double PW92_dalpha_drs(double rs);

// GGA / Meta GGA helpers
double e_X_ueg(double rho);
double eps_c_pw92(double rho_a, double rho_b);
double ke_density_ueg(double rho);

// VV10
double VV10_kernel(double b, double C, double R2, double rho_1, double rho_2, double nrm_grho_1, double nrm_grho_2);

// Old
double density2(double x, double y, double z, const Molecule& mol, const Matrix& P);

#endif
