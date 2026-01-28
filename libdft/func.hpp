#ifndef FUNCHEADERDEF
#define FUNCHEADERDEF

#include "../libint/1e.hpp"
#include "dft_helper.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>

class XC{
	public:
		XC(const std::string& method);
		
		void (*xc_functional)(XC*);
		bool restricted;

		Matrix* P   = nullptr;
		Matrix* FXC = nullptr;		// restricted

		Matrix* P_A   = nullptr;	// unrestricted
		Matrix* P_B   = nullptr;	//
		Matrix* FXC_A = nullptr;	//
		Matrix* FXC_B = nullptr;	//

		const std::vector<std::vector<std::vector<std::vector<double>>>>* eris = nullptr;
		const Molecule* mol = nullptr;
		const grid* g = nullptr;

		double E_XC = 0.0;
		
		double rho   = 0.0;
		double rho_a = 0.0;
		double rho_b = 0.0;

		std::vector<double> gradient_rho   = {0.0, 0.0, 0.0};
		std::vector<double> gradient_rho_a = {0.0, 0.0, 0.0};
		std::vector<double> gradient_rho_b = {0.0, 0.0, 0.0};

		double ke_density   = 0.0;
		double ke_density_a = 0.0;
		double ke_density_b = 0.0;
};

// *_ret objects are per-grid point properties

struct LDA_ret{
	double e_XC;				// energy density
	std::vector<double> v_XC;	// potential(s, a/b)
};

struct GGA_ret{
	double e_XC;									// 	energy density
	std::vector<double> drho_XC;					// 	\frac{de_{XC}^{GGA}}{d\rho_{\sigma}}
	std::vector<double> dgamma_XC; 					/*	\frac{de_{XC}^{GGA}}{d\gamma_{\sigma\sigma'}}

	gamma_{\sigma\sigma'} = \nabla \rho_{\sigma} \cdot \nabla \rho_{\sigma'}
	if restricted: dgamma_XC[0] is derivative w.r.t. gamma = |\nabla \rho|^2
	else: dgamma_XC is derivative w.r.t. {\gamma_{aa}, \gamma_{bb}, \gamma_{ab}} */
};

struct MGGA_ret{
	double e_XC;
	std::vector<double> drho_XC;
	std::vector<double> dgamma_XC;
	std::vector<double> dtau_XC;
};

extern std::unordered_map<std::string, void (*)(XC*)> xc_register;
void call_xc_functional(XC* xc);

// HF //
void R_HF_X(XC* xc);
void U_HF_X(XC* xc);
void R_HF_SNX(XC* xc);	// Seminumerical Exchange

// LDA //
void Slater(XC* xc);
void VWN5(XC* xc);
void PW92(XC* xc);

// GGA //
void PBE_X(XC* xc);
void PBE(XC* xc);
/*
// Meta GGA //
template <int nlcflag, typename F>
XC_ret F_XC_MGGA(const XC_inp& inp, F&& v_MGGA){
	Molecule* m = inp.mol;
	grid*     g = inp.g;
	Matrix*  pa = inp.PA;
	Matrix*  pb = inp.PB;
	
	double rho_a, rho_b;
	std::vector<double> grho_a(3);
	std::vector<double> grho_b(3);
	double t_a, t_b;
	
	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(m->AOs.size());
	std::vector<double> gpx_buf(m->AOs.size());
	std::vector<double> gpy_buf(m->AOs.size());
	std::vector<double> gpz_buf(m->AOs.size());
	std::vector<double> tmp_grd(3);
	
	Matrix F_XC_A(pa->rows, pa->cols);
	Matrix F_XC_B(pb->rows, pb->cols);

	for(int i = 0; i < num_gpts; i++){
		for(int j = 0; j < m->AOs.size(); j++){
			phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			tmp_grd = m->AOs[j].evaluate_gradient(g->x[i], g->y[i], g->z[i]);
			gpx_buf[j] = tmp_grd[0];
			gpy_buf[j] = tmp_grd[1];
			gpz_buf[j] = tmp_grd[2];
		}
		rho_a = density(phi_buf, *pa);
		rho_b = density(phi_buf, *pb);
		grho_a = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
		grho_b = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
		t_a = ke_density(gpx_buf, gpy_buf, gpz_buf, *pa);
		t_b = ke_density(gpx_buf, gpy_buf, gpz_buf, *pb);
		
		for(int mu = 0; mu < F_XC_A.rows; mu++){
			F_XC_A.matrix[mu][mu] += g->w[i] * v_MGGA(rho_a, rho_b, grho_a, grho_b, t_a, t_b, 0);
			F_XC_B.matrix[mu][mu] += g->w[i] * v_MGGA(rho_a, rho_b, grho_a, grho_b, t_a, t_b, 1);
			for(int nu = 0; nu < mu; nu++){
				F_XC_A.matrix[mu][nu] += g->w[i] * v_MGGA(rho_a, rho_b, grho_a, grho_b, t_a, t_b, 0);
				F_XC_A.matrix[nu][mu] = F_XC_A.matrix[mu][nu];
				F_XC_B.matrix[mu][nu] += g->w[i] * v_MGGA(rho_a, rho_b, grho_a, grho_b, t_a, t_b, 1);
				F_XC_B.matrix[nu][mu] = F_XC_B.matrix[mu][nu];
			}
		}
	}
	return {F_XC_A, F_XC_B};
}

// nlcflag = 0 : no nonlocal correlation; = 1 : nonlocal correlation (like VV10...)
template <int nlcflag, typename F>
double E_XC_MGGA(const XC_inp& inp, F&& e_MGGA){
	Molecule* m = inp.mol;
	grid*     g = inp.g;
	Matrix*  pa = inp.PA;
	Matrix*  pb = inp.PB;

	double rho_a, rho_b;
	std::vector<double> grho_a(3);
	std::vector<double> grho_b(3);
	double tau_a, tau_b;
	
	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(m->AOs.size());
	std::vector<double> gpx_buf(m->AOs.size());
	std::vector<double> gpy_buf(m->AOs.size());
	std::vector<double> gpz_buf(m->AOs.size());
	std::vector<double> tmp_grd(3);

	double E_XC = 0;
	for(int i = 0; i < num_gpts; i++){
		for(int j = 0; j < m->AOs.size(); j++){
			phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			tmp_grd = m->AOs[j].evaluate_gradient(g->x[i], g->y[i], g->z[i]);
			gpx_buf[j] = tmp_grd[0];
			gpy_buf[j] = tmp_grd[1];
			gpz_buf[j] = tmp_grd[2];
		}
		rho_a = density(phi_buf, *pa);
		rho_b = density(phi_buf, *pb);
		grho_a = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
		grho_b = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
		tau_a = ke_density(gpx_buf, gpy_buf, gpz_buf, *pa);
		tau_b = ke_density(gpx_buf, gpy_buf, gpz_buf, *pb);

		if constexpr (nlcflag == 0){
			E_XC += g->w[i] * e_MGGA(rho_a, rho_b, grho_a, grho_b, tau_a, tau_b);
		}
		else if constexpr (nlcflag == 1){
			E_XC += g->w[i] * e_MGGA(rho_a, rho_b, grho_a, grho_b, tau_a, tau_b, m, g, pa, pb, i);
		}
		else {assert((nlcflag == 0) || (nlcflag == 1));}
	}
	return E_XC;
}

XC_ret U_B97M_V(const XC_inp& inp);
double U_B97M_V_E(const XC_inp& inp);
*/
#endif
