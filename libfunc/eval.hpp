#ifndef EVALHEADERDEF
#define EVALHEADERDEF

#include "../libmath/linalg.hpp"

#include <omp.h>
#include <vector>

class  XC;
struct XC_inp;
struct XC_ret;

void scf_xc_call(XC* xc);

// Evaluate E_XC, F_XC
void LDA (XC* xc);
void GGA_MGGA(XC* xc);

// Helpers
void zero_xc_data(XC* inp, int spins);

void eval_bfs_per_gpt(const XC& xc, std::vector<double>& phi_buf, int gpix);
void eval_bfs_per_gpt(const XC& xc, Matrix& phi_buf, int gpix);
void eval_bfs_per_batch(const XC& xc, Matrix& phi_buf, int g_start, int g_end);
void eval_bfs_grad_per_gpt(const XC& xc, std::vector<double>& phi_buf, std::vector<double>& gpx_buf, std::vector<double>& gpy_buf, 
	std::vector<double>& gpz_buf, std::vector<double>& temp_grad, int gpix);

void eval_density_per_gpt(const XC& xc, XC_inp& inp, const std::vector<double>& phi_buf);
void eval_density_grad_per_gpt(const XC& xc, XC_inp& inp, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf);
void eval_density_grad_ke_per_gpt(const XC& xc, XC_inp& inp, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf);

double LDA_per_gpt(const XC& xc, const XC_inp& inp, XC_ret& ret, const std::vector<double>& phi_buf, int spins, int gpix);
double GGA_per_gpt(const XC& xc, const XC_inp& inp, XC_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int spins, int gpix);
double MGGA_per_gpt(const XC& xc, const XC_inp& inp, XC_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int spins, int gpix);

double GGA_F_second_term(const XC_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, const std::vector<double>* grho[2], int mu, int nu, int s);
double MGGA_F_third_term(const XC_ret& ret, const std::vector<double>& gpx_buf, const std::vector<double>& gpy_buf, 
	const std::vector<double>& gpz_buf, int mu, int nu, int s);

#endif
