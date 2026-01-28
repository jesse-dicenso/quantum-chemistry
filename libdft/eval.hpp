#ifndef EVALHEADERDEF
#define EVALHEADERDEF

#include "../libmath/linalg.hpp"

#include <vector>

class  XC;
struct LDA_ret;
struct GGA_ret;
struct MGGA_ret;

// Evaluate E_XC, F_XC
void LDA (XC* xc, LDA_ret  (*func)(XC*));
void GGA (XC* xc, GGA_ret  (*func)(XC*));
void MGGA(XC* xc, MGGA_ret (*func)(XC*));

// Helpers
void zero_xc_data(XC* inp);

void eval_bfs_per_gpt(XC* xc, std::vector<double>& phi_buf, int gpix);
void eval_bfs_per_gpt(XC* xc, Matrix& phi_buf, int gpix);
void eval_bfs_grad_per_gpt(XC* xc, std::vector<double>& phi_buf, std::vector<double>& gpx_buf, std::vector<double>& gpy_buf, 
	std::vector<double>& gpz_buf, std::vector<double>& temp_grad, int gpix);

void eval_density_per_gpt(XC* xc, const std::vector<double>& phi_buf);
void eval_density_grad_per_gpt(XC* xc, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf);
void eval_density_grad_ke_per_gpt(XC* xc, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf);

void LDA_per_gpt(XC* xc, LDA_ret (*func)(XC*), const std::vector<double>& phi_buf, int gpix);
void GGA_per_gpt(XC* xc, GGA_ret (*func)(XC*), const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int gpix);
void MGGA_per_gpt(XC* xc, MGGA_ret (*func)(XC*), const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int gpix);

double GGA_F_second_term(const GGA_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, std::vector<double>* grho[2], int mu, int nu, int s);
double MGGA_F_third_term(const MGGA_ret& ret, const std::vector<double>& gpx_buf, const std::vector<double>& gpy_buf, 
	const std::vector<double>& gpz_buf, int mu, int nu, int s);

#endif
