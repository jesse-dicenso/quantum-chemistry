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

		int main_iter = 0;

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

// Meta GGA //
void B97M_V(XC* xc);

#endif
