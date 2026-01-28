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
	{ "R_HF_SNX" , R_HF_SNX},
	{ "R_Slater" , Slater  },
	{ "U_Slater" , Slater  },
	{ "R_VWN5"   , VWN5    },
	{ "U_VWN5"   , VWN5    },
	{ "R_PW92"   , PW92    },
	{ "U_PW92"   , PW92    },
	{ "R_PBE_X"  , PBE_X   },
	{ "U_PBE_X"  , PBE_X   },
	{ "R_PBE"    , PBE     },
	{ "U_B97M-V" , B97M_V  }
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

// Helper for *_HF_SNX
void SNX_A(XC* xc, Matrix& A, int gpix){
	const std::vector<GF>& bfs = xc->mol->AOs;
	const std::vector<double> xyzg = {xc->g->x[gpix], xc->g->y[gpix], xc->g->z[gpix]};
	for(int sg = 0; sg < A.rows; sg++){
		for(int nu = 0; nu < A.cols; nu++){
			A.matrix[sg][nu] = V(bfs[sg], bfs[nu], xyzg);
		}
	}
}

void R_HF_SNX(XC* xc){
	assert((xc->P!=nullptr) && (xc->g!=nullptr) && (xc->mol!=nullptr) && (xc->FXC!=nullptr));
	const int K = xc->mol->AOs.size();
	const int gpts = xc->g->num_gridpoints;
	const std::vector<double>& w = xc->g->w;
	const Matrix* p = xc->P;
	Matrix *fxc = xc->FXC;
	Matrix X(K, 1);
	Matrix A(K, K);
	Matrix G(K, 1);
	zero_xc_data(xc);
	for(int g = 0; g < gpts; g++){
		eval_bfs_per_gpt(xc, X, g);
		SNX_A(xc, A, g);
		G = (A * (*p * X)) * w[g];
		for(int mu = 0; mu < K; mu++){
			for(int nu = 0; nu < K; nu++){
				fxc->matrix[mu][nu] -= 0.5 * X.matrix[mu][0] * G.matrix[nu][0];
				xc->E_XC += (g==(gpts-1) ? fxc->matrix[mu][nu] * p->matrix[mu][nu] : 0.0);
			}
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
	inline const double RX = -cbrt(3.0 / M_PI);
	inline const double UX = -cbrt(6.0 / M_PI);
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
		ret.e_XC = RX * (3.0 / 4.0) * rho * rho_3;
		ret.v_XC = {RX * rho_3};
		return ret;
	};
	auto func_u = [](XC* inp) {
		const double rho_a   = inp->rho_a;
		const double rho_b   = inp->rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);	
		LDA_ret ret;
		ret.e_XC = UX * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		ret.v_XC = { UX * rho_a_3, UX * rho_b_3 };
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
		const double e_X   = RX * (3.0 / 4.0) * rho * rho_3;
		const double v_X   = RX * rho_3;

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
		const double e_X     = UX * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		std::vector<double> v = {UX * rho_a_3, UX * rho_b_3};
		
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
		const double e_X   = RX * (3.0 / 4.0) * rho * rho_3;
		const double v_X   = RX * rho_3;
		
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
		const double e_X     = UX * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		std::vector<double> v = {UX * rho_a_3, UX * rho_b_3};
		
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
		const double e_LDA = RX * (3.0 / 4.0) * rho * rho_3;
		const double v_LDA = RX * rho_3;

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
		const double e_LDA_a = (3.0 / 4.0) * RX * rho_a * rho_a_3;
		const double e_LDA_b = (3.0 / 4.0) * RX * rho_b * rho_b_3;
		const double v_LDA_a = RX * rho_a_3;	// v_LDA_X spin scaling already accounted for!
		const double v_LDA_b = RX * rho_b_3;	// v_LDA_X spin scaling already accounted for!
		
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
	auto func = [](XC* inp) {
		// Slater Exchange
		const double rho = inp->rho;
		const double rho_3 = cbrt(rho);
		const double e_LDA_X = RX * (3.0 / 4.0) * rho * rho_3;
		const double v_LDA_X = RX * rho_3;

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
	GGA(xc, func);
}

// MGGA ///////////////////////////////////////////////////////

namespace _B97M_V{
	inline const double c_tau_ueg = 3.0 * cbrt(36 * M_PI) * M_PI / 5.0;

	inline constexpr double c_x00 = 1.000;
	inline constexpr double c_x10 = 0.416;
	inline constexpr double c_x01 = 1.308;
	inline constexpr double c_x11 = 3.070;
	inline constexpr double c_x02 = 1.901;

	inline constexpr double c_css00 =  1.000;
	inline constexpr double c_css10 = -5.668;
	inline constexpr double c_css02 = -1.855;
	inline constexpr double c_css32 = -20.497;
	inline constexpr double c_css42 = -20.364;

	inline constexpr double c_cos00 =  1.000;
	inline constexpr double c_cos10 =  2.535;
	inline constexpr double c_cos01 =  1.573;
	inline constexpr double c_cos32 = -6.427;
	inline constexpr double c_cos03 = -6.298;

	inline constexpr double gamma_x   = 0.004;
	inline constexpr double gamma_css = 0.2;
	inline constexpr double gamma_cos = 0.006;

	inline constexpr double b_VV10 = 6.00;
	inline constexpr double C_VV10 = 0.01;
}

void B97M_V_E(XC* xc){
	using namespace _B97M_V;
	assert(!xc->restricted);
	assert((xc->mol!=nullptr) && (xc->g!=nullptr) && (xc->P_A!=nullptr) && (xc->P_B!=nullptr));	
	/* Separation of same/opp spin correlation 
	   requires unrestricted calculation!
	   E = e_X + e_css + e_cos + e_VV10 */
	auto func = [](XC* inp)
	{
		const double rho   = inp->rho;
		const double rho_a = inp->rho_a;
		const double rho_b = inp->rho_b;
		const std::vector<double>& grho   = &inp->gradient_rho;
		const std::vector<double>& grho_a = &inp->gradient_rho_a;
		const std::vector<double>& grho_b = &inp->gradient_rho_b;
		const double tau_a = inp->ke_density_a;
		const double tau_b = inp->ke_density_b;

		const double nrm_grho   = sqrt(grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2]);
		const double nrm_grho_a = sqrt(grho_a[0] * grho_a[0] + grho_a[1] * grho_a[1] + grho_a[2] * grho_a[2]);
		const double nrm_grho_b = sqrt(grho_b[0] * grho_b[0] + grho_b[1] * grho_b[1] + grho_b[2] * grho_b[2]);

		const double tau_ueg_a = c_tau_ueg * rho_a * cbrt(rho_a * rho_a);
		const double tau_ueg_b = c_tau_ueg * rho_b * cbrt(rho_b * rho_b);

		const double t_a = tau_ueg_a / tau_a;
		const double t_b = tau_ueg_b / tau_b;
		const double dt_a_drho_a = (5.0 / 3.0) * t_a / rho_a;
		const double dt_b_drho_b = (5.0 / 3.0) * t_b / rho_b;
		const double dt_a_dtau_a = -t_a / tau_a;
		const double dt_b_dtau_b = -t_b / tau_b;

		// e_X
		const double e_X_ueg_a = UX * (3.0 / 4.0) * (rho_a * cbrt(rho_a));
		const double e_X_ueg_b = UX * (3.0 / 4.0) * (rho_b * cbrt(rho_b));
		const double s_a = nrm_grho_a / cbrt(intpow(rho_a, 4));
		const double s_b = nrm_grho_b / cbrt(intpow(rho_b, 4));
		const double s2_a = s_a * s_a;
		const double s2_b = s_b * s_b;
		const double wx_a = (t_a - 1) / (t_a + 1);
		const double wx_b = (t_b - 1) / (t_b + 1);
		const double ux_a = gamma_x * s2_a / (1 + gamma_x * s2_a);
		const double ux_b = gamma_x * s2_b / (1 + gamma_x * s2_b);

		const double gx_a = c_x00 + (c_x10 * wx_a) + (c_x01 * ux_a) + (c_x11 * wx_a * ux_a) + (c_x02 * ux_a * ux_a);
		const double gx_b = c_x00 + (c_x10 * wx_b) + (c_x01 * ux_b) + (c_x11 * wx_b * ux_b) + (c_x02 * ux_b * ux_b);

		const double e_X = e_X_ueg_a * gx_a + e_X_ueg_b * gx_b;

		// e_X derivatives
		const double de_X_ueg_a_drho_a = UX * cbrt(rho_a);
		const double de_X_ueg_b_drho_b = UX * cbrt(rho_b);
		const double ds2_a_drho_a = -(8.0 / 3.0) * s2_a / rho_a;
		const double ds2_b_drho_b = -(8.0 / 3.0) * s2_b / rho_b;
		const double ds2_a_dgrho2_a = s2_a / (nrm_grho_a * nrm_grho_a);
		const double ds2_b_dgrho2_b = s2_b / (nrm_grho_b * nrm_grho_b);
		const double dwx_a_dt_a = 2.0 / ((t_a + 1) * (t_a + 1));
		const double dwx_b_dt_b = 2.0 / ((t_b + 1) * (t_b + 1));
		const double dux_a_ds2_a = gamma_x / ((1 + gamma_x * s2_a) * (1 + gamma_x * s2_a));
		const double dux_b_ds2_b = gamma_x / ((1 + gamma_x * s2_b) * (1 + gamma_x * s2_b));
		const double dgx_a_dwx_a = c_x10 + (c_x11 * ux_a);	
		const double dgx_b_dwx_b = c_x10 + (c_x11 * ux_b);
		const double dgx_a_dux_a = c_x01 + (c_x11 * wx_a) + (2 * c_x02 * ux_a);
		const double dgx_b_dux_b = c_x01 + (c_x11 * wx_b) + (2 * c_x02 * ux_b);
	
		const double de_X_drho_a = de_X_ueg_drho_a * gx_a + 
			e_X_ueg_a * (dgx_a_dwx_a * dwx_a_dt_a * dt_a_drho_a + dgx_a_dux_a * dux_a_ds2_a * ds2_a_drho_a);
		const double de_X_drho_b = de_X_ueg_drho_b * gx_b + 
			e_X_ueg_b * (dgx_b_dwx_b * dwx_b_dt_b * dt_b_drho_b + dgx_b_dux_b * dux_b_ds2_b * ds2_b_drho_b);

		const double de_X_dgrho2_a = e_X_ueg_a * dgx_a_dux_a * dux_a_ds2_a * ds2_a_dgrho2_a;
		const double de_X_dgrho2_b = e_X_ueg_b * dgx_a_dux_b * dux_a_ds2_b * ds2_b_dgrho2_b;

		const double de_X_dtau_a = e_X_ueg_a * dgx_a_dwx_a * dwx_a_dt_a * dt_a_dtau_a;
		const double de_X_dtau_b = e_X_ueg_b * dgx_b_dwx_b * dwx_b_dt_b * dt_b_dtau_b;
	
		// e_css
		const double eps_pw92_aa = eps_c_pw92(rho_a, 0.0);
		const double eps_pw92_bb = eps_c_pw92(0.0, rho_b);
		const double wc_aa = wx_a;
		const double wc_bb = wx_b;
		const double uc_aa = gamma_css * s_a * s_a / (1 + gamma_css * s_a * s_a);
		const double uc_bb = gamma_css * s_b * s_b / (1 + gamma_css * s_b * s_b);

		const double gcss_aa = c_css00 + (c_css10 * wc_aa) + (c_css02 * uc_aa * uc_aa) + (c_css32 * intpow(wc_aa, 3) * uc_aa * uc_aa) + 
			(c_css42 * intpow(wc_aa, 4) * uc_aa * uc_aa);
		const double gcss_bb = c_css00 + (c_css10 * wc_bb) + (c_css02 * uc_bb * uc_bb) + (c_css32 * intpow(wc_bb, 3) * uc_bb * uc_bb) + 
			(c_css42 * intpow(wc_bb, 4) * uc_bb * uc_bb);

		const double e_css = (rho_a * eps_pw92_aa * gcss_aa) + (rho_b * eps_pw92_bb * gcss_bb);

		// e_css_derivatives
		const double deps_pw92_aa_drho_a = deps_c_dns_pw92(rho_a, 0.0, 0);
		const double deps_pw92_bb_drho_b = deps_c_dns_pw92(0.0, rho_b, 1);
		const double dwc_aa_dt_a = dwx_a_dt_a;
		const double dwc_bb_dt_b = dwx_b_dt_b;
		const double duc_ab_ds2_a = gamma_css / ((1 + gamma_css * s2_a) * (1 + gamma_css * s2_a));
		const double duc_bb_ds2_b = gamma_css / ((1 + gamma_css * s2_b) * (1 + gamma_css * s2_b));
		
		const double dgcss_aa_dwc_aa = c_css10 + (3 * c_css32 * wc_aa * wc_aa * uc_aa * uc_aa) + 
			(4 * c_css42 * intpow(wc_aa, 3) * uc_aa * uc_aa);
		const double dgcss_bb_dwc_bb = c_css10 + (3 * c_css32 * wc_bb * wc_bb * uc_bb * uc_bb) + 
			(4 * c_css42 * intpow(wc_bb, 3) * uc_bb * uc_bb);

		const double dgcss_aa_duc_aa = (2 * c_css02 * uc_aa) + (2 * c_css32 * intpow(wc_aa, 3) * uc_aa) + 
			(2 * c_css42 * intpow(wc_aa, 4) * uc_aa); 
		const double dgcss_bb_duc_bb = (2 * c_css02 * uc_bb) + (2 * c_css32 * intpow(wc_bb, 3) * uc_bb) + 
			(2 * c_css42 * intpow(wc_bb, 4) * uc_bb);

		const double de_css_drho_a = (eps_pw92_aa + rho_a * deps_pw92_drho_a) * gcss_aa + 
			rho_a * eps_pw92_aa * (dgcss_aa_dwc_aa * dwc_aa_dt_a * dt_a_drho_a + dgcss_aa_duc_aa * duc_aa_ds2_a * ds2_a_drho_a);
		const double de_css_drho_b = (eps_pw92_bb + rho_b * deps_pw92_drho_b) * gcss_bb + 
			rho_b * eps_pw92_bb * (dgcss_bb_dwc_bb * dwc_bb_dt_b * dt_b_drho_b + dgcss_bb_duc_bb * duc_bb_ds2_b * ds2_b_drho_b);
	
		const double de_css_dgrho2_a = rho_a * eps_pw92_aa * (dgcss_aa_duc_aa * duc_aa_ds2_a * ds2_a_dgrho2_a);
		const double de_css_dgrho2_b = rho_b * eps_pw92_bb * (dgcss_bb_duc_bb * duc_bb_ds2_b * ds2_b_dgrho2_b);

		const double de_css_dtau_a = rho_a * eps_pw92_aa * (dgcss_aa_dwc_aa * dwc_aa_dt_a * dt_a_dtau_a);
		const double de_css_dtau_b = rho_b * eps_pw92_bb * (dgcss_bb_dwc_bb * dwc_bb_dt_b * dt_b_dtau_b);
	
		// e_cos
		const double t_ab = 0.5 * (t_a + t_b);
		const double s2ab = 0.5 * (s_a * s_a + s_b * s_b);
		const double e_pw92_ab = rho * eps_c_pw92(rho_a, rho_b) - rho_a * eps_pw92_aa - rho_b * eps_pw92_bb;
		const double wc_ab = (t_ab - 1) / (t_ab + 1);
		const double uc_ab = gamma_cos * s2ab / (1 + gamma_cos * s2ab);

		const double gcos = c_cos00 + (c_cos10 * wc_ab) + (c_cos01 * uc_ab) + (c_cos32 * intpow(wc_ab, 3) * uc_ab * uc_ab) + 
			(c_cos03 * intpow(uc_ab, 3));

		// e_cos derivatives
		const double de_pw92_ab_drho_a = (eps_c_pw92(rho_a, rho_b) + rho * deps_c_dns_pw92(rho_a, rho_b, 0)) - 
			(eps_pw92_aa + rho_a * deps_pw92_aa_drho_a);
		const double de_pw92_ab_drho_b = (eps_c_pw92(rho_a, rho_b) + rho * deps_c_dns_pw92(rho_a, rho_b, 1)) - 
			(eps_pw92_bb + rho_b * deps_pw92_bb_drho_b);
		const double dwc_ab_dt_ab = 2.0 / ((t_ab + 1) * (t_ab + 1));
		const double duc_ab_ds2_ab = gamma_cos / ((1 + gamma_cos * s2ab) * (1 + gamma_cos * s2ab));
		
		const double dgcos_dwc_ab = c_cos10 + (3 * c_cos32 * wc_ab * wc_ab * uc_ab * uc_ab);
		const double dgcos_duc_ab = c_cos01 + (2 * c_cos32 * wc_ab * wc_ab * wc_ab * uc_ab) + (3 * c_cos03 * uc_ab * uc_ab);

		const double de_cos_drho_a = de_pw92_ab_drho_a * gcos + e_pw92_ab * 
			(dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_a_drho_a + dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_a_drho_a);
		const double de_cos_drho_b = de_pw92_ab_drho_b * gcos + e_pw92_ab * 
			(dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_b_drho_b + dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_b_drho_b);

		const double de_cos_dgrho2_a = e_pw92_ab * (dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_a_dgrho2_a);
		const double de_cos_dgrho2_b = e_pw92_ab * (dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_b_dgrho2_b);

		const double de_cos_dtau_a = e_pw92_ab * (dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_a_dtau_a);
		const double de_cos_dtau_b = e_pw92_ab * (dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_b_dtau_b);

		// VV10

		GGA_ret VV10_ret = VV10_per_gpt(*xc, rho, grho2, b_VV10, C_VV10);
		const double e_VV10 = VV10_ret.e_XC;
		const double de_VV10_drho = VV10_ret.drho_XC[0];
		const double de_VV10_dgrho2 = VV10_ret.dgamma_XC[0];

		MGGA_ret ret;
		ret.e_XC = e_X + e_css + e_cos + e_VV10;
		ret.drho_XC = {
			(de_X_drho_a + de_css_drho_a + de_cos_drho_a + de_VV10_drho),
			(de_X_drho_b + de_css_drho_b + de_cos_drho_b + de_VV10_drho)
		};
		ret.dgamma_XC = {
			(de_X_dgrho2_a + de_css_dgrho2_a + de_cos_dgrho2_a + de_VV10_dgrho2), 
			(de_X_dgrho2_b + de_css_dgrho2_b + de_cos_dgrho2_b + de_VV10_dgrho2)
		};
		ret.dtau_XC = {
			(de_X_dtau_a + de_css_dtau_a + de_cos_dtau_a), 
			(de_X_dtau_b + de_css_dtau_b + de_cos_dtau_b)
		};
		return ret;
	};
	MGGA(xc, func);
}

GGA_ret VV10_per_gpt(const XC& xc, double ref_rho, double ref_grho2, const double b, const double C){
	using namespace _B97M_V;
	
	const int ref_gpt = xc.main_iter;
	const double beta = sqrt(sqrt(27.0) / b) / 32.0 / b;
	
	double e_VV10  = 0.0;
	double F_rho   = 0.0;
	double F_gamma = 0.0;
	std::vector<double> phi_buf(m->AOs.size());
	std::vector<double> gpx_buf(m->AOs.size());
	std::vector<double> gpy_buf(m->AOs.size());
	std::vector<double> gpz_buf(m->AOs.size());
	std::vector<double> tmp_grd(3);
	double rho_gpt, grho2_gpt, R2, PHI;

	Matrix* p = xc->P_T;
	const std::vector<double>& gx = xc.g->x;
	const std::vector<double>& gy = xc.g->y;
	const std::vector<double>& gz = xc.g->z;
	const std::vector<double>& gw = xc.g->w;
	for(int gpt = 0; gpt < xc.g->num_gridpoints; gpt++){
		eval_bfs_grad_per_gpt(&xc, phi_buf, gpx_buf, gpy_buf, gpz_buf, tmp_grd, gpt);
		rho_gpt = density(phi_buf, *p);
		tmp_grd = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *p);
		grho2_gpt = (
			tmp_grd[0] * tmp_grd[0] +
			tmp_grd[1] * tmp_grd[1] +
			tmp_grd[2] * tmp_grd[2]
		);	
		R2 = intpow(x[ref_gpt] - x[gpt], 2) + intpow(y[ref_gpt] - y[gpt], 2) + intpow(z[ref_gpt] - z[gpt], 2);
		PHI = /* VV10_kernel */;
		e_VV10 += gw[gpt] * rho_gpt * PHI;
	}
	e_VV10 = ref_rho * (beta + 0.5 * e_VV10);
	F_rho += beta;

	GGA_ret ret;
	ret.e_XC = e_VV10;
	ret.drho_XC = {F_rho};
	ret.dgamma_XC = {F_gamma};
	return ret;
}

