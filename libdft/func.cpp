#include "func.hpp"
#include "eval.hpp"

XC::XC(const std::string& method){
	if		(method.substr(0,2)=="R_"){restricted=true;}
	else if (method.substr(0,2)=="U_"){restricted=false;}
	else{throw std::invalid_argument("ERR: method must be restricted 'R_' or unrestricted 'U_'");}
	xc_functional = xc_register[method];
}

void call_xc_functional(XC* xc){
	xc->xc_functional(xc);
}

std::unordered_map<std::string, void (*)(XC*)> xc_register = 
{
	{ "R_HF"     , R_HF_X  },
	{ "U_HF"     , U_HF_X  },
	{ "R_Slater" , Slater  },
	{ "U_Slater" , Slater  },
	{ "R_VWN5"   , VWN5    },
	{ "U_VWN5"   , VWN5    },
	{ "R_PW92"   , PW92    },
	{ "U_PW92"   , PW92    },
	{ "R_PBE_X"  , PBE_X   },
	{ "U_PBE_X"  , PBE_X   },
	{ "R_PBE"    , PBE     }/*,
	{ "U_B97M-V" , B97M_V  },*/
};

// HF //
void R_HF_X(XC* xc){
	assert((xc->P!=nullptr) && (xc->eris!=nullptr) && (xc->FXC!=nullptr));
	double fxc = 0;
	xc->E_XC = 0;
	for(int mu = 0; mu < xc->FXC->rows; mu++){
		for(int nu = 0; nu < xc->FXC->cols; nu++){
			fxc = 0;
			for(int ld = 0; ld < xc->FXC->rows; ld++){
				for(int sg = 0; sg < xc->FXC->cols; sg++){
					fxc -= xc->P->matrix[ld][sg] * (*xc->eris)[mu][ld][sg][nu];
				}
			}
			xc->FXC->matrix[mu][nu] = 0.5 * fxc;
			xc->E_XC += xc->FXC->matrix[mu][nu] * xc->P->matrix[mu][nu];
		}
	}
	xc->E_XC *= 0.5;
}

void U_HF_X(XC* xc){
	assert((xc->P_A!=nullptr) && (xc->P_B!=nullptr) && (xc->eris!=nullptr) && (xc->FXC_A!=nullptr) && (xc->FXC_B!=nullptr));
	double fxc_a, fxc_b;
	xc->E_XC = 0;
	for(int mu = 0; mu < xc->FXC_A->rows; mu++){
		for(int nu = 0; nu < xc->FXC_A->cols; nu++){
			fxc_a = 0;
			fxc_b = 0;
			for(int ld = 0; ld < xc->FXC_A->rows; ld++){
				for(int sg = 0; sg < xc->FXC_A->cols; sg++){
					fxc_a -= xc->P_A->matrix[ld][sg] * (*xc->eris)[mu][ld][sg][nu];
					fxc_b -= xc->P_B->matrix[ld][sg] * (*xc->eris)[mu][ld][sg][nu];
				}
			}
			xc->FXC_A->matrix[mu][nu] = fxc_a;
			xc->FXC_B->matrix[mu][nu] = fxc_b;
			xc->E_XC += xc->P_A->matrix[mu][nu] * xc->FXC_A->matrix[mu][nu] + xc->P_B->matrix[mu][nu] * xc->FXC_B->matrix[mu][nu];
		}
	}
	xc->E_XC *= 0.5;
}

///////////////////////////////////////////////////////////////
// !!!													 !!! //
// !!!   Lambdas below give PER GRID-POINT F_XC / E_XC   !!! //
// !!!													 !!! //
///////////////////////////////////////////////////////////////

// LDA ////////////////////////////////////////////////////////

namespace _SLATER{
	inline const double R = -cbrt(3.0 / M_PI);
	inline const double U = -cbrt(6.0 / M_PI);
}

void Slater(XC* xc){
	using namespace _SLATER;
	assert((xc->mol!=nullptr) && (xc->g!=nullptr));
	if(xc->restricted){assert(xc->P!=nullptr);}
	else{assert((xc->P_A!=nullptr) && (xc->P_B!=nullptr));}

	auto func_r = [](XC* inp) {
		const double rho   = inp->rho;
		const double rho_3 = cbrt(rho);
		LDA_ret ret;
		ret.e_XC = R * (3.0 / 4.0) * rho * rho_3;
		ret.v_XC = {R * rho_3};
		return ret;
	};
	auto func_u = [](XC* inp) {
		const double rho_a   = inp->rho_a;
		const double rho_b   = inp->rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);	
		LDA_ret ret;
		ret.e_XC = U * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		ret.v_XC = { U * rho_a_3, U * rho_b_3 };
		return ret;
	};
	LDA_ret (*func)(XC*) = (xc->restricted ? func_r : func_u);
	LDA(xc, func);
}

namespace _VWN5{
	// Paramagnetic (zeta = 0)
	inline const     double A_P  = (1 - log(2)) / (M_PI * M_PI);
	inline constexpr double x0_P = -0.10498;
	inline constexpr double b_P  =  3.72744;
	inline constexpr double c_P  =  12.9352;
	inline constexpr double X0_P =  x0_P * x0_P + b_P * x0_P + c_P;
	inline const     double Q_P  =  sqrt(4 * c_P - b_P * b_P);
	// Ferromagnetic (zeta = 1)
	inline const     double A_F  =  A_P / 2;
	inline constexpr double x0_F = -0.32500;
	inline constexpr double b_F  =  7.06042;
	inline constexpr double c_F  =  18.0578;
	inline constexpr double X0_F =  x0_F * x0_F + b_F * x0_F + c_F;
	inline const     double Q_F  =  sqrt(4 * c_F - b_F * b_F);
	inline const     double ddf0 = 4.0 / (9.0 * (cbrt(2) - 1));
}

void VWN5(XC* xc){
	using namespace _SLATER;
	using namespace _VWN5;
	assert((xc->mol!=nullptr) && (xc->g!=nullptr));
	if(xc->restricted){assert(xc->P!=nullptr);}
	else{assert((xc->P_A!=nullptr) && (xc->P_B!=nullptr));}
	
	auto func_r = [](XC* inp) {
		const double rho = inp->rho;
		const double rho_3 = cbrt(rho);
		const double e_X   = R * (3.0 / 4.0) * rho * rho_3;
		const double v_X   = R * rho_3;

		const double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		const double X  = x * x + b_P * x + c_P;

		const double eps_c = A_P * (
			log(x * x / X) + (2 * b_P / Q_P) * (1 - (2 * x0_P + b_P) * x0_P / X0_P) * atan(Q_P / (2 * x + b_P)) - 
			(b_P * x0_P / X0_P) * log((x - x0_P) * (x - x0_P) / X)
		);
		const double v_c = eps_c - A_P * (x / (3 * X)) * (c_P / x - b_P * x0_P / (x - x0_P));

		LDA_ret ret;
		ret.e_XC = e_X + rho * eps_c;
		ret.v_XC = { v_X + v_c };
		return ret;
	};
	auto func_u = [](XC* inp) {
		const double rho_a   = inp->rho_a;
		const double rho_b   = inp->rho_b;
		const double rho     = rho_a + rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);
		const double e_X     = U * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		std::vector<double> v = {U * rho_a_3, U * rho_b_3};
		
		const double x = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		const double zeta = ( rho_a - rho_b ) / rho;
		const double zeta3 = zeta * zeta * zeta;
		const std::vector<double> dzeta_drho = {2 * rho_b / (rho * rho), -2 * rho_a / (rho * rho)};
		const double f  = f_zeta(zeta);
		const double df = df_zeta(zeta);
		const double alpha = VWN_alpha(x);
		const double dalpha_drho = VWN_dalpha_drho(x, rho);
		const double X_P  = x * x + b_P * x + c_P;
		const double X_F  = x * x + b_F * x + c_F;
		
		const double eps_c_P = A_P * (
			log(x * x / X_P) + (2 * b_P / Q_P) * (1 - (2 * x0_P + b_P) * x0_P / X0_P) * 
			atan(Q_P / (2 * x + b_P)) - (b_P * x0_P / X0_P) * log((x - x0_P) * (x - x0_P) / X_P)
		);
		const double eps_c_F = A_F * (
			log(x * x / X_F) + (2 * b_F / Q_F) * (1 - (2 * x0_F + b_F) * x0_F / X0_F) * 
			atan(Q_F / (2 * x + b_F)) - (b_F * x0_F / X0_F) * log((x - x0_F) * (x - x0_F) / X_F)
		);
		const double e_c = rho * eps_c_P + rho * alpha * (f / ddf0) * (1 - zeta3 * zeta) + rho * (eps_c_F - eps_c_P) * f * zeta3 * zeta;
		
		const double v_c_P = eps_c_P - A_P * (x / (3 * X_P)) * (c_P / x - b_P * x0_P / (x - x0_P));
		const double v_c_F = eps_c_F - A_F * (x / (3 * X_F)) * (c_F / x - b_F * x0_F / (x - x0_F));

		for(int s = 0; s < 2; s++){
			v[s] += v_c_P + (alpha + rho * dalpha_drho) * (f / ddf0) * (1 - zeta3 * zeta) + 
					rho * alpha * ((df/ddf0) * (1 - zeta3 * zeta) - 4 * zeta3 * (f / ddf0)) * dzeta_drho[s] +
					(v_c_F - v_c_P) * f * zeta3 * zeta + 
					rho * (eps_c_F - eps_c_P) * (df * zeta3 * zeta + 4 * zeta3 * f) * dzeta_drho[s];	
		}
		LDA_ret ret;
		ret.e_XC = e_X + e_c;
		ret.v_XC = v;
		return ret;
	};
	LDA_ret (*func)(XC*) = (xc->restricted ? func_r : func_u);
	LDA(xc, func);
}

namespace _PW92 {
	// Paramagnetic (zeta = 0)
	inline const     double A_P  = (1 - log(2)) / (M_PI * M_PI);
	inline constexpr double a1_P = 0.21370;
	inline constexpr double b1_P = 7.5957;
	inline constexpr double b2_P = 3.5876;
	inline constexpr double b3_P = 1.6382;
	inline constexpr double b4_P = 0.49294;
	// Ferromagnetic (zeta = 1)
	inline const     double A_F  = A_P / 2;
	inline constexpr double a1_F = 0.20548;
	inline constexpr double b1_F = 14.1189;
	inline constexpr double b2_F = 6.1977;
	inline constexpr double b3_F = 3.3662;
	inline constexpr double b4_F = 0.62517;
}

void PW92(XC* xc){
	using namespace _SLATER;
	using namespace _PW92;
	assert((xc->mol!=nullptr) && (xc->g!=nullptr));
	if(xc->restricted){assert(xc->P!=nullptr);}
	else{assert((xc->P_A!=nullptr) && (xc->P_B!=nullptr));}

	auto func_r = [](XC* inp) {
		const double rho = inp->rho;
		const double rho_3 = cbrt(rho);
		const double e_X   = R * (3.0 / 4.0) * rho * rho_3;
		const double v_X   = R * rho_3;
		
		const double rs = cbrt(3 / (4 * M_PI * rho));
		const double Q0  = -2 * A_P * (1 + a1_P * rs);
		const double Q1  =  2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs);
		const double Q1p =      A_P * (b1_P / sqrt(rs) + 2 * b2_P + 3 * b3_P * sqrt(rs) + 4 * b4_P * rs);
		
		const double eps_c = -2 * A_P * (1 + a1_P * rs) * log(1 + 
			1 / (2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs)));
		const double v_c = eps_c - (rs / 3) * (-2 * A_P * a1_P * log(1 + 1 / Q1) - Q0 * Q1p / (Q1 * Q1 + Q1));
		
		LDA_ret ret;
		ret.e_XC = e_X + rho * eps_c;
		ret.v_XC = { v_X + v_c };
		return ret;
	};
	auto func_u = [](XC* inp) {
		const double rho_a   = inp->rho_a;
		const double rho_b   = inp->rho_b;
		const double rho     = rho_a + rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);
		const double e_X     = U * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		std::vector<double> v = {U * rho_a_3, U * rho_b_3};
		
		double rs = cbrt(3 / (4 * M_PI * rho));
		const double zeta = (rho_a - rho_b) / rho;
		const double zeta3 = zeta * zeta * zeta;
		const double zeta4 = zeta3 * zeta;
		const double f = f_zeta(zeta);
		const double df = df_zeta(zeta);
		const double alpha = PW92_alpha(rs);
		const double dalpha_drs = PW92_dalpha_drs(rs);
		
		const double eps_P = -2 * A_P * (1 + a1_P * rs) * log(1 + 1 / (2 * A_P * 
			(b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs)));
		const double eps_F = -2 * A_F * (1 + a1_F * rs) * log(1 + 1 / (2 * A_F * 
			(b1_F * sqrt(rs) + b2_F * rs + b3_F * sqrt(intpow(rs, 3)) + b4_F * rs * rs)));
		const double eps = eps_P + alpha * (f / _VWN5::ddf0) * (1 - zeta4) + (eps_F - eps_P) * f * zeta4;
		
		const double Q0_P   = -2 * A_P * (1 + a1_P * rs);
		const double Q1_P   =  2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs);
		const double Q1p_P  =      A_P * (b1_P / sqrt(rs) + 2 * b2_P + 3 * b3_P * sqrt(rs) + 4 * b4_P * rs);
		const double deps_P = -2 * A_P * a1_P * log(1 + 1 / Q1_P) - Q0_P * Q1p_P / (Q1_P * Q1_P + Q1_P);
	
		const double Q0_F   = -2 * A_F * (1 + a1_F * rs);
		const double Q1_F   =  2 * A_F * (b1_F * sqrt(rs) + b2_F * rs + b3_F * sqrt(intpow(rs, 3)) + b4_F * rs * rs);
		const double Q1p_F  =      A_F * (b1_F / sqrt(rs) + 2 * b2_F + 3 * b3_F * sqrt(rs) + 4 * b4_F * rs);
		const double deps_F = -2 * A_F * a1_F * log(1 + 1 / Q1_F) - Q0_F * Q1p_F / (Q1_F * Q1_F + Q1_F);
		
		const double deps_dr = deps_P * (1 - f * zeta4) + deps_F * f * zeta4 + dalpha_drs * (f / _VWN5::ddf0) * (1 - zeta4);
		const double deps_dz = 4 * zeta3 * f * (eps_F - eps_P - alpha / _VWN5::ddf0) + 
			df * (zeta4 * (eps_F - eps_P) + (1 - zeta4) * alpha / _VWN5::ddf0);
		
		v[0] += eps - (rs / 3) * deps_dr - (zeta - 1) * deps_dz;
		v[1] += eps - (rs / 3) * deps_dr - (zeta + 1) * deps_dz;

		LDA_ret ret;
		ret.e_XC = e_X + rho * eps;
		ret.v_XC = v;
		return ret;
	};
	LDA_ret (*func)(XC*) = (xc->restricted ? func_r : func_u);
	LDA(xc, func);
}

// GGA ////////////////////////////////////////////////////////

namespace _PBE{
	inline constexpr double beta = 0.066725;
	inline constexpr double kappa = 0.804;
	inline constexpr double mu = beta * (M_PI * M_PI / 3);
	inline const     double gamma = (1 - log(2)) / (M_PI * M_PI);
}

void PBE_X(XC* xc){
	using namespace _SLATER;
	using namespace _PBE;
	assert((xc->mol!=nullptr) && (xc->g!=nullptr));
	if(xc->restricted){assert(xc->P!=nullptr);}
	else{assert((xc->P_A!=nullptr) && (xc->P_B!=nullptr));}
	
	auto func_r = [](XC* inp) {
		// Slater Exchange
		const double rho = inp->rho;
		const double rho_3 = cbrt(rho);
		const double e_LDA = R * (3.0 / 4.0) * rho * rho_3;
		const double v_LDA = R * rho_3;

		// Enhancement Factor
		const double grho2 = (
			inp->gradient_rho[0] * inp->gradient_rho[0] + 
			inp->gradient_rho[1] * inp->gradient_rho[1] + 
			inp->gradient_rho[2] * inp->gradient_rho[2] 
		);
		const double kF = cbrt(3 * M_PI * M_PI * rho);
		const double s2 = grho2 / (4 * kF * kF * rho * rho);
		const double ds2_drho = -8.0 * s2 / (3.0 * rho);
		const double ds2_dgrho2 = s2 / grho2;	
		const double FX_d = 1 + mu * s2 / kappa;	
		const double FX = 1 + kappa - kappa / FX_d;
		const double dFX_ds2 = mu / (FX_d * FX_d);
		const double dFX_drho = dFX_ds2 * ds2_drho; 
		const double dFX_dgrho2 = dFX_ds2 * ds2_dgrho2;
 
		const double de_drho = v_LDA * FX + e_LDA * dFX_drho;
		const double de_dgrho2 = e_LDA * dFX_dgrho2;

		GGA_ret ret;
		ret.e_XC = e_LDA * FX;
		ret.drho_XC = {de_drho};
		ret.dgamma_XC = {de_dgrho2};
		return ret;
	};
	auto func_u = [](XC* inp) {
		// Slater Exchange
		const double rho   = 2 * inp->rho;
		const double rho_a = 2 * inp->rho_a;
		const double rho_b = 2 * inp->rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);
		const double e_LDA_a = (3.0 / 4.0) * R * rho_a * rho_a_3;
		const double e_LDA_b = (3.0 / 4.0) * R * rho_b * rho_b_3;
		const double v_LDA_a = R * rho_a_3;	// v_LDA_X spin scaling already accounted for!
		const double v_LDA_b = R * rho_b_3;	// v_LDA_X spin scaling already accounted for!
		
		// Enhancement Factor
		const double grho2_a = 4 * (
			inp->gradient_rho_a[0] * inp->gradient_rho_a[0] + 
			inp->gradient_rho_a[1] * inp->gradient_rho_a[1] + 
			inp->gradient_rho_a[2] * inp->gradient_rho_a[2]
		); 
		const double grho2_b = 4 * (
			inp->gradient_rho_b[0] * inp->gradient_rho_b[0] + 
			inp->gradient_rho_b[1] * inp->gradient_rho_b[1] + 
			inp->gradient_rho_b[2] * inp->gradient_rho_b[2]
		); 
		const double kF_a = cbrt(3 * M_PI * M_PI * rho_a);
		const double kF_b = cbrt(3 * M_PI * M_PI * rho_b);
		const double s2_a = grho2_a / (4 * kF_a * kF_a * rho_a * rho_a);
		const double s2_b = grho2_b / (4 * kF_b * kF_b * rho_b * rho_b);
		const double ds2_drho_a = -8.0 * s2_a / (3.0 * rho_a);
		const double ds2_drho_b = -8.0 * s2_b / (3.0 * rho_b);
		const double ds2_dgrho2_a = s2_a / grho2_a;	
		const double ds2_dgrho2_b = s2_b / grho2_b;	
		const double FX_d_a = 1 + mu * s2_a / kappa;	
		const double FX_d_b = 1 + mu * s2_b / kappa;	
		const double FX_a = 1 + kappa - kappa / FX_d_a;
		const double FX_b = 1 + kappa - kappa / FX_d_b;
		const double dFX_ds2_a = mu / (FX_d_a * FX_d_a);
		const double dFX_ds2_b = mu / (FX_d_b * FX_d_b);
		const double dFX_drho_a = dFX_ds2_a * ds2_drho_a;
		const double dFX_drho_b = dFX_ds2_b * ds2_drho_b;
		const double dFX_dgrho2_a = dFX_ds2_a * ds2_dgrho2_a;
		const double dFX_dgrho2_b = dFX_ds2_b * ds2_dgrho2_b;
		
		const double de_drho_a = v_LDA_a * FX_a + e_LDA_a * dFX_drho_a;	// All divided by
		const double de_drho_b = v_LDA_b * FX_b + e_LDA_b * dFX_drho_b;	// 2 to account
		const double de_dgrho2_a = 2 * e_LDA_a * dFX_dgrho2_a;			// for spin scaling
		const double de_dgrho2_b = 2 * e_LDA_b * dFX_dgrho2_b;			// of e_X and v_X

		GGA_ret ret;
		ret.e_XC = 0.5 * (e_LDA_a * FX_a + e_LDA_b * FX_b);
		ret.drho_XC = {de_drho_a, de_drho_b};
		ret.dgamma_XC = {de_dgrho2_a, de_dgrho2_b};
		return ret;
	};
	GGA_ret (*func)(XC*) = (xc->restricted ? func_r : func_u);
	GGA(xc, func);
}

void PBE(XC* xc){
	using namespace _SLATER;
	using namespace _PW92;
	using namespace _PBE;
	assert(xc->restricted); // only restricted PBE for now!
	assert((xc->mol!=nullptr) && (xc->g!=nullptr) && (xc->P!=nullptr));	
	auto func_r = [](XC* inp) {
		// Slater Exchange
		const double rho = inp->rho;
		const double rho_3 = cbrt(rho);
		const double e_LDA_X = R * (3.0 / 4.0) * rho * rho_3;
		const double v_LDA_X = R * rho_3;

		// Enhancement Factor
		const double grho2 = (
			inp->gradient_rho[0] * inp->gradient_rho[0] + 
			inp->gradient_rho[1] * inp->gradient_rho[1] + 
			inp->gradient_rho[2] * inp->gradient_rho[2] 
		);
		const double kF = cbrt(3 * M_PI * M_PI * rho);
		const double s2 = grho2 / (4 * kF * kF * rho * rho);
		const double ds2_drho = -8.0 * s2 / (3.0 * rho);
		const double ds2_dgrho2 = s2 / grho2;	
		const double FX_d = 1 + mu * s2 / kappa;	
		const double FX = 1 + kappa - kappa / FX_d;
		const double dFX_ds2 = mu / (FX_d * FX_d);
		const double dFX_drho = dFX_ds2 * ds2_drho; 
		const double dFX_dgrho2 = dFX_ds2 * ds2_dgrho2;
 
		const double de_X_drho = v_LDA_X * FX + e_LDA_X * dFX_drho;
		const double de_X_dgrho2 = e_LDA_X * dFX_dgrho2;

		GGA_ret ret;
		ret.e_XC = e_LDA_X * FX;
		ret.drho_XC = {de_X_drho};
		ret.dgamma_XC = {de_X_dgrho2};

		// PW92 correlation
		const double rs = cbrt(3 / (4 * M_PI * rho));
		const double Q0  = -2 * A_P * (1 + a1_P * rs);
		const double Q1  =  2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs);
		const double Q1p =      A_P * (b1_P / sqrt(rs) + 2 * b2_P + 3 * b3_P * sqrt(rs) + 4 * b4_P * rs);
		
		const double eps_c_LDA = -2 * A_P * (1 + a1_P * rs) * log(1 + 
			1 / (2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs)));
		const double deps_c_LDA_dn = -(rs / 3) * (-2 * A_P * a1_P * log(1 + 1 / Q1) - Q0 * Q1p / (Q1 * Q1 + Q1)) / rho;
		const double v_c_LDA = eps_c_LDA + rho * deps_c_LDA_dn;

		// PBE correlation correction
		const double ks = sqrt(4 * kF / M_PI);
		const double t2 = grho2 / (4 * ks * ks * rho * rho);
		double A_PBE = (beta / gamma) / (exp(-eps_c_LDA / gamma) - 1);
		A_PBE = (A_PBE > 1e10 ? 1e10 : A_PBE);
		const double dnm = 1 + A_PBE * t2 + A_PBE * A_PBE * t2 * t2;
		const double Q = 1 + (beta / gamma) * t2 * (1 + A_PBE * t2) / dnm;
		const double H = gamma * log(Q);

		const double dH_dt2 = (beta / Q) * (1 + 2 * A_PBE * t2) / (dnm * dnm);
		const double dt2_dn = -(7.0 / 3.0) * t2 / rho;
		const double dH_dA_PBE = -(beta / Q) * A_PBE * t2 * t2 * t2 * (2 + A_PBE * t2) / (dnm * dnm);
		const double dA_PBE_deps = A_PBE * (A_PBE + beta / gamma) / beta;
		const double dt2_dgrho2 = t2 / grho2;

		ret.e_XC += rho * (eps_c_LDA + H);
		ret.drho_XC[0] += v_c_LDA + H + rho * (dH_dt2 * dt2_dn + dH_dA_PBE * dA_PBE_deps * deps_c_LDA_dn);
		ret.dgamma_XC[0] += rho * dH_dt2 * dt2_dgrho2;
		return ret;
	};
	GGA_ret (*func)(XC*) = func_r;
	GGA(xc, func);
}

/*
// MGGA ///////////////////////////////////////////////////////
Matrix F_VV10(double rho_a, double rho_b, const std::vector<double>& grho_a, const std::vector<double>& grho_b, 
			  double phi1, double phi2, double gpx1, double gpy1, double gpz1, double gpx2, double gpy2, double gpz2, 
			  Molecule* m, grid* g, Matrix* pa, Matrix* pb, int gidx)
{
	const double b = 6.00;
	const double C = 0.01;

	Matrix F(pa->rows, pa->cols);
	double rho = rho_a + rho_b;
	if (rho < 1e-16) {return 0.0;}
	double nrm_grho = sqrt(grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2]);

	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(m->AOs.size());
	std::vector<double> gpx_buf(m->AOs.size());
	std::vector<double> gpy_buf(m->AOs.size());
	std::vector<double> gpz_buf(m->AOs.size());
	std::vector<double> tmp_grd(3);
	double rho_a_j, rho_b_j, rho_j, nrm_grho_j, R2;
	std::vector<double> grho_a_j(3);
	std::vector<double> grho_b_j(3);
	std::vector<double> grho_j(3);
	for(int j = 0; j < num_gpts; j++){
		for(int k = 0; k < m->AOs.size(); k++){
			phi_buf[k] = m->AOs[k].evaluate(g->x[j], g->y[j], g->z[j]);
			tmp_grd = m->AOs[k].evaluate_gradient(g->x[j], g->y[j], g->z[j]);
			gpx_buf[k] = tmp_grd[0];
			gpy_buf[k] = tmp_grd[1];
			gpz_buf[k] = tmp_grd[2];
		}

		R2 = intpow(g->x[gidx] - g->x[j], 2) + intpow(g->y[gidx] - g->y[j], 2) + intpow(g->z[gidx] - g->z[j], 2);

		rho_a_j = density(phi_buf, *pa);
		rho_b_j = density(phi_buf, *pb);
		rho_j = rho_a_j + rho_b_j;
		grho_a_j = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
		grho_b_j = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
		
		grho_j[0] = grho_a_j[0] + grho_b_j[0];
		grho_j[1] = grho_a_j[1] + grho_b_j[1];
		grho_j[2] = grho_a_j[2] + grho_b_j[2];

		nrm_grho_j = sqrt(grho_j[0] * grho_j[0] + grho_j[1] * grho_j[1] + grho_j[2] * grho_j[2]);

		e_VV10 += g->w[j] * VV10_kernel(b, C, R2, rho, rho_j, nrm_grho, nrm_grho_j) * rho_j;
	}
	e_VV10 *= 0.5;
	e_VV10 += sqrt(sqrt(intpow(3 / (b * b), 3))) / 32;
	e_VV10 *= rho;
}

XC_ret U_B97M_V(const XC_inp& inp){
	// Parameters
	double c_x00, c_x10, c_x01, c_x11, c_x02;
	double c_css00, c_css10, c_css02, c_css32, c_css42;
	double c_cos00, c_cos10, c_cos01, c_cos32, c_cos03;
	double b, C;

	c_x00 = 1.000;
	c_x10 = 0.416;
	c_x01 = 1.308;
	c_x11 = 3.070;
	c_x02 = 1.901;

	c_css00 =  1.000;
	c_css10 = -5.668;
	c_css02 = -1.855;
	c_css32 = -20.497;
	c_css42 = -20.364;

	c_cos00 =  1.000;
	c_cos10 =  2.535;
	c_cos01 =  1.573;
	c_cos32 = -6.427;
	c_cos03 = -6.298;

	b = 6.00;
	C = 0.01;

	auto v_ = [c_x00, c_x10, c_x01, c_x11, c_x02, 
			  c_css00, c_css10, c_css02, c_css32, c_css42, 
			  c_cos00, c_cos10, c_cos01, c_cos32, c_cos03, 
			  ]
	(double rho_a, double rho_b, const std::vector<double>& grho_a, const std::vector<double>& grho_b, double tau_a, double tau_b, int spin)
	{
		//
	};

	return F_XC_MGGA<1>(inp, v);
}

double U_B97M_V_E(const XC_inp& inp){
	// Parameters
	double c_x00, c_x10, c_x01, c_x11, c_x02;
	double c_css00, c_css10, c_css02, c_css32, c_css42;
	double c_cos00, c_cos10, c_cos01, c_cos32, c_cos03;
	double b, C;

	c_x00 = 1.000;
	c_x10 = 0.416;
	c_x01 = 1.308;
	c_x11 = 3.070;
	c_x02 = 1.901;

	c_css00 =  1.000;
	c_css10 = -5.668;
	c_css02 = -1.855;
	c_css32 = -20.497;
	c_css42 = -20.364;

	c_cos00 =  1.000;
	c_cos10 =  2.535;
	c_cos01 =  1.573;
	c_cos32 = -6.427;
	c_cos03 = -6.298;

	b = 6.00;
	C = 0.01;
	
	// E = e_X + e_css + e_cos + e_VV10
	auto e = [c_x00, c_x10, c_x01, c_x11, c_x02, 
			  c_css00, c_css10, c_css02, c_css32, c_css42, 
			  c_cos00, c_cos10, c_cos01, c_cos32, c_cos03, 
			  b, C]
	(double rho_a, double rho_b, const std::vector<double>& grho_a, const std::vector<double>& grho_b, double tau_a, double tau_b, 
	 Molecule* m, grid* g, Matrix* pa, Matrix* pb, int gidx)
	{
		double rho = rho_a + rho_b;
		if(rho < 1e-20){return 0.0;}
		std::vector<double> grho = {grho_a[0] + grho_b[0], grho_a[1] + grho_b[1], grho_a[2] + grho_b[2]};

		double nrm_grho_a = sqrt(grho_a[0] * grho_a[0] + grho_a[1] * grho_a[1] + grho_a[2] * grho_a[2]);
		double nrm_grho_b = sqrt(grho_b[0] * grho_b[0] + grho_b[1] * grho_b[1] + grho_b[2] * grho_b[2]);
		double nrm_grho   = sqrt(grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2]);

		double t_a = ke_density_ueg(rho_a) / tau_a;
		double t_b = ke_density_ueg(rho_b) / tau_b;
		
		// e_X
		double s_a, s_b, wx_a, wx_b, ux_a, ux_b, gx_a, gx_b, e_X;
		s_a = nrm_grho_a / cbrt(intpow(rho_a, 4));
		s_b = nrm_grho_b / cbrt(intpow(rho_b, 4));
		wx_a = (t_a - 1) / (t_a + 1);
		wx_b = (t_b - 1) / (t_b + 1);
		ux_a = 0.004 * s_a * s_a / (1 + 0.004 * s_a * s_a);
		ux_b = 0.004 * s_b * s_b / (1 + 0.004 * s_b * s_b);

		gx_a = c_x00 + c_x10 * wx_a + c_x01 * ux_a + c_x11 * wx_a * ux_a + c_x02 * ux_a * ux_a;
		gx_b = c_x00 + c_x10 * wx_b + c_x01 * ux_b + c_x11 * wx_b * ux_b + c_x02 * ux_b * ux_b;

		e_X = e_X_ueg(rho_a) * gx_a + e_X_ueg(rho_b) * gx_b;
		
		// e_css
		double e_pw92_aa, e_pw92_bb, wc_aa, wc_bb, uc_aa, uc_bb, gcss_aa, gcss_bb, e_css;
		e_pw92_aa = rho_a * eps_c_pw92(rho_a, 0);
		e_pw92_bb = rho_b * eps_c_pw92(0, rho_b);
		wc_aa = wx_a;
		wc_bb = wx_b;
		uc_aa = 0.2 * s_a * s_a / (1 + 0.2 * s_a * s_a);
		uc_bb = 0.2 * s_b * s_b / (1 + 0.2 * s_b * s_b);

		gcss_aa = c_css00 + c_css10 * wc_aa + c_css02 * uc_aa * uc_aa + c_css32 * intpow(wc_aa, 3) * uc_aa * uc_aa + 
			 	  c_css42 * intpow(wc_aa, 4) * uc_aa * uc_aa;
		gcss_bb = c_css00 + c_css10 * wc_bb + c_css02 * uc_bb * uc_bb + c_css32 * intpow(wc_bb, 3) * uc_bb * uc_bb + 
				  c_css42 * intpow(wc_bb, 4) * uc_bb * uc_bb;

		e_css = e_pw92_aa * gcss_aa + e_pw92_bb * gcss_bb;
		
		// e_cos
		double t_ab, s2ab, e_pw92_ab, wc_ab, uc_ab, gcos, e_cos;
		t_ab = 0.5 * (t_a + t_b);
		s2ab = 0.5 * (s_a * s_a + s_b * s_b);
		e_pw92_ab = rho * eps_c_pw92(rho_a, rho_b) - e_pw92_aa - e_pw92_bb;
		wc_ab = (t_ab - 1) / (t_ab + 1);
		uc_ab = 0.006 * s2ab / (1 + 0.006 * s2ab);

		gcos = c_cos00 + c_cos10 * wc_ab + c_cos01 * uc_ab + c_cos32 * intpow(wc_ab, 3) * uc_ab * uc_ab + c_cos03 * intpow(uc_ab, 3);

		e_cos = e_pw92_ab * gcos;

		// e_VV10	
		double e_VV10 = 0;

		int num_gpts = g->num_gridpoints;
		std::vector<double> phi_buf(m->AOs.size());
		std::vector<double> gpx_buf(m->AOs.size());
		std::vector<double> gpy_buf(m->AOs.size());
		std::vector<double> gpz_buf(m->AOs.size());
		std::vector<double> tmp_grd(3);
		double rho_a_j, rho_b_j, rho_j, nrm_grho_j, R2;
		std::vector<double> grho_a_j(3);
		std::vector<double> grho_b_j(3);
		std::vector<double> grho_j(3);
		for(int j = 0; j < num_gpts; j++){
			for(int k = 0; k < m->AOs.size(); k++){
				phi_buf[k] = m->AOs[k].evaluate(g->x[j], g->y[j], g->z[j]);
				tmp_grd = m->AOs[k].evaluate_gradient(g->x[j], g->y[j], g->z[j]);
				gpx_buf[k] = tmp_grd[0];
				gpy_buf[k] = tmp_grd[1];
				gpz_buf[k] = tmp_grd[2];
			}

			R2 = intpow(g->x[gidx] - g->x[j], 2) + intpow(g->y[gidx] - g->y[j], 2) + intpow(g->z[gidx] - g->z[j], 2);
	
			rho_a_j = density(phi_buf, *pa);
			rho_b_j = density(phi_buf, *pb);
			rho_j = rho_a_j + rho_b_j;
			grho_a_j = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
			grho_b_j = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
			
			grho_j[0] = grho_a_j[0] + grho_b_j[0];
			grho_j[1] = grho_a_j[1] + grho_b_j[1];
			grho_j[2] = grho_a_j[2] + grho_b_j[2];

			nrm_grho_j = sqrt(grho_j[0] * grho_j[0] + grho_j[1] * grho_j[1] + grho_j[2] * grho_j[2]);

			e_VV10 += g->w[j] * VV10_kernel(b, C, R2, rho, rho_j, nrm_grho, nrm_grho_j) * rho_j;
		}
		e_VV10 *= 0.5;
		e_VV10 += sqrt(sqrt(intpow(3 / (b * b), 3))) / 32;
		e_VV10 *= rho;

		return e_X + e_css + e_cos + e_VV10;
	};
	return E_XC_MGGA<1>(inp, e);
}
*/
