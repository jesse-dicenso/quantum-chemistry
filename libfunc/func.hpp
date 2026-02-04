#ifndef FUNCHEADERDEF
#define FUNCHEADERDEF

#include "../libint/1e.hpp"
#include "dft_helper.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>

// XC_* contains per-grid point quantities
struct XC_inp{
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

struct XC_ret{
	double e_XC;                    // energy density
	std::vector<double> drho_XC;    // \frac{de_{XC}}{d\rho_{\sigma}}
	std::vector<double> dgamma_XC;  // \frac{de_{XC}}{d\gamma_{\sigma\sigma'}}
	std::vector<double> dtau_XC;    // \frac{de_{XC}}{d\tau_{\sigma}}
    // restricted   : arr[0] only
    // unrestricted : arr[0] -> aa, arr[1] -> bb, arr[2] -> ab if mixed terms required
    std::vector<Matrix> F_XC;       // {fxc} if restricted, {a, b} if unrestricted
};

class XC{
	public:
		XC(const std::string& method);
		
		void (*xc_functional)(const XC& xc, const XC_inp&, XC_ret& ret) = nullptr;
        void (*nlc_functional)(XC* xc) = nullptr;
		bool restricted;
        bool isHF   = false;
        bool isHFSN = false;
        bool isLDA  = false;
        bool isGGA  = false;
        bool isMGGA = false;
        bool isNLC  = false;

        std::vector<Matrix*> P;
        std::vector<Matrix*> F_XC;

		const std::vector<std::vector<std::vector<std::vector<double>>>>* eris = nullptr;
		const Molecule* mol = nullptr;
		const grid* g = nullptr;

        double E_XC = 0.0;

        // For VV10
        const Matrix* overlap = nullptr;
        std::vector<double> nlc_params;
};

extern std::unordered_map<std::string, unsigned int> xc_register; 

// HF //
void HFX  (XC* xc);
void HFSNX(XC* xc);	// Seminumerical Exchange

// LDA //
void Slater(const XC& xc, const XC_inp& inp, XC_ret& ret);
void VWN5  (const XC& xc, const XC_inp& inp, XC_ret& ret);
void PW92  (const XC& xc, const XC_inp& inp, XC_ret& ret);

// GGA //
void PBE_X(const XC& xc, const XC_inp& inp, XC_ret& ret);
void PBE  (const XC& xc, const XC_inp& inp, XC_ret& ret);

// Meta GGA //
void B97M_V(const XC& xc, const XC_inp& inp, XC_ret& ret);

// NLC //
void VV10(XC* xc);
//XC_ret VV10_per_gpt(XC* xc, double ref_rho, double ref_grho2, const double b, const double C);

#endif
