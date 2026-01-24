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

std::unordered_map<std::string, void(XC*)> xc_register = 
{
	{ "R_HF"     , R_HF_X   },
	{ "U_HF"     , U_HF_X   },
	{ "R_Slater" , Slater },
	{ "U_Slater" , Slater } /*,
	{ "R_VWN5"   , R_VWN5   },
	{ "U_VWN5"   , U_VWN5   },
	{ "R_PW92"   , R_PW92   },
	{ "U_PW92"   , U_PW92   },
	{ "R_PBE"    , R_PBE    },
	{ "U_PBE_X"  , U_PBE_X  },
	{ "U_B97M-V" , B97M_V   },*/
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
// !!! Functionals below give PER GRID-POINT F_XC / E_XC !!! //
// !!!													 !!! //
///////////////////////////////////////////////////////////////

// LDA ////////////////////////////////////////////////////////

void Slater(XC* xc){
	assert((xc.mol!=nullptr) && (xc.g!=nullptr));
	if(xc->restricted){assert(xc->P!=nullptr);}
	else{assert((xc->P_A!=nullptr) && (xc->P_B!=nullptr));}

	constexpr C_X = -(3.0 / 4.0) * cbrt(3.0 / M_PI);
	auto func_r = [C_X](XC* inp) {
		XC_ret ret;
		double rho_3 = cbrt(xc->rho);

		ret.E_XC =  C_X * intpow(rho_3, 4);
		ret.V_XC = {C_X * (4.0 / 3.0) * rho_3};

		return xc_ret;
	};
	auto func_u = [C_X](XC* inp) {
		XC_ret ret;
		double rho_a_3 = cbrt(xc->rho_a);
		double rho_b_3 = cbrt(xc->rho_b);
		
		ret.E_XC =  C_X * cbrt(2.0) * (intpow(rho_a_3, 4) + intpow(rho_b_3, 4));
		ret.V_XC = {C_X * cbrt(2.0) * (4.0 / 3.0) * rho_a_3, 
					C_X * cbrt(2.0) * (4.0 / 3.0) * rho_b_3 };
 
		return xc_ret;
	};
	std::vector<double>(XC*) func = (xc->restricted ? func_r : func_u);
	LDA(xc, func);
}
/*

double R_Slater(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	auto e = [](double rho) {
		return cbrt(rho * rho * rho * rho);
	};
	return -(3.0/4.0) * cbrt(3.0 / M_PI) * E_XC_LDA<0>(inp, e);
}

XC_ret U_Slater_X(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	auto v = [](double rho_s) {
		return -cbrt(6 * rho_s / M_PI);
	};
	return F_XC_LDA<1>(inp, v);
}

double U_Slater_X_E(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	auto e = [](double rho_a, double rho_b) {
		return cbrt(rho_a * rho_a * rho_a * rho_a) + cbrt(rho_b * rho_b * rho_b * rho_b);
	};
	return -(3.0/4.0) * cbrt(6.0 / M_PI) * E_XC_LDA<1>(inp, e);
}

XC_ret R_VWN5_c(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double x0 = -0.10498;
	const double b  =  3.72744;
	const double c  =  12.9352;
	const double X0 = x0 * x0 + b * x0 + c;
	const double Q  = sqrt(4 * c - b * b);
	auto v = [A, x0, b, c, X0, Q](double rho) {
		if(rho < 1e-20) {return 0.0;}
		double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		double X  = x * x + b * x + c;
		double vc = (
			log(x * x / X) + (2 * b / Q) * (1 - (2 * x0 + b) * x0 / X0) * atan(Q / (2 * x + b)) - (b * x0 / X0) * log((x - x0) * (x - x0) / X)
		);
		vc -= (x / (3 * X)) * (c / x - b * x0 / (x - x0));
		vc *= A;
		return vc; 
	};
	return F_XC_LDA<0>(inp, v);
}

double R_VWN5_c_E(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double x0 = -0.10498;
	const double b  =  3.72744;
	const double c  =  12.9352;
	const double Q  = sqrt(4 * c - b * b);
	const double X0 = x0 * x0 + b * x0 + c;
	auto e = [A, x0, b, c, X0, Q](double rho) {
		if(rho < 1e-20) {return 0.0;}
		double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		double X  = x * x + b * x + c;
		double ec = rho * A * (
			log(x * x / X) + (2 * b / Q) * (1 - (2 * x0 + b) * x0 / X0) * atan(Q / (2 * x + b)) - (b * x0 / X0) * log((x - x0) * (x - x0) / X)
		);
		return ec; 
	};
	return E_XC_LDA<0>(inp, e);
}

XC_ret U_VWN5_c(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	// zeta = 0 constants
	const double A_0  = (1 - log(2)) / (M_PI * M_PI);
	const double x0_0 = -0.10498;
	const double b_0  =  3.72744;
	const double c_0  =  12.9352;
	const double X0_0 =  x0_0 * x0_0 + b_0 * x0_0 + c_0;
	const double Q_0  =  sqrt(4 * c_0 - b_0 * b_0);
	// zeta = 1 constants
	const double A_1  =  A_0 / 2;
	const double x0_1 = -0.32500;
	const double b_1  =  7.06042;
	const double c_1  =  18.0578;
	const double X0_1 =  x0_1 * x0_1 + b_1 * x0_1 + c_1;
	const double Q_1  =  sqrt(4 * c_1 - b_1 * b_1);
	auto v = [A_0, x0_0, b_0, c_0, X0_0, Q_0, A_1, x0_1, b_1, c_1, X0_1, Q_1](double rho_a, double rho_b, int spin) 
	{
		double rho = rho_a + rho_b;
		if(rho < 1e-20) {return 0.0;}
		double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));

		double zeta = ( rho_a - rho_b ) / rho;
		double zeta3 = zeta * zeta * zeta;
		// spin = 0 -> alpha, spin = 1 -> beta
		double dzeta_drho;
		if(spin == 0) {dzeta_drho = 2 * rho_b / (rho * rho);}
		else if(spin == 1) {dzeta_drho = -2 * rho_a / (rho * rho);}
		else{assert((spin==0) || (spin==1));}
		double f  = f_zeta(zeta);
		double df = df_zeta(zeta);
		double ddf0 = 4.0 / ( 9.0 * ( cbrt(2) - 1 ) );
		double alpha = VWN_alpha(x);
		double dalpha_drho = VWN_dalpha_drho(x, rho);

		// evaluate energy densities and potentials for zeta = 0, 1
		double X_0  = x * x + b_0 * x + c_0;
		double X_1  = x * x + b_1 * x + c_1;

		double ec_0 = rho * A_0 * (
			log(x * x / X_0) + (2 * b_0 / Q_0) * (1 - (2 * x0_0 + b_0) * x0_0 / X0_0) * 
			atan(Q_0 / (2 * x + b_0)) - (b_0 * x0_0 / X0_0) * log((x - x0_0) * (x - x0_0) / X_0)
		);
		double ec_1 = rho * A_1 * (
			log(x * x / X_1) + (2 * b_1 / Q_1) * (1 - (2 * x0_1 + b_1) * x0_1 / X0_1) * 
			atan(Q_1 / (2 * x + b_1)) - (b_1 * x0_1 / X0_1) * log((x - x0_1) * (x - x0_1) / X_1)
		);
	
		double vc_0 = (
			log(x * x / X_0) + (2 * b_0 / Q_0) * (1 - (2 * x0_0 + b_0) * x0_0 / X0_0) * 
			atan(Q_0 / (2 * x + b_0)) - (b_0 * x0_0 / X0_0) * log((x - x0_0) * (x - x0_0) / X_0)
		);
		vc_0 -= (x / (3 * X_0)) * (c_0 / x - b_0 * x0_0 / (x - x0_0));
		vc_0 *= A_0;	
		double vc_1 = (
			log(x * x / X_1) + (2 * b_1 / Q_1) * (1 - (2 * x0_1 + b_1) * x0_1 / X0_1) * 
			atan(Q_1 / (2 * x + b_1)) - (b_1 * x0_1 / X0_1) * log((x - x0_1) * (x - x0_1) / X_1)
		);
		vc_1 -= (x / (3 * X_1)) * (c_1 / x - b_1 * x0_1 / (x - x0_1));
		vc_1 *= A_1;

		double vc_s = vc_0 + (alpha + rho * dalpha_drho) * (f / ddf0) * (1 - zeta3 * zeta) + 
					  rho * alpha * ((df/ddf0) * (1 - zeta3 * zeta) - 4 * zeta3 * (f / ddf0)) * dzeta_drho +
					  (vc_1 - vc_0) * f * zeta3 * zeta + 
					  (ec_1 - ec_0) * (df * zeta3 * zeta + 4 * zeta3 * f) * dzeta_drho;
					  
		return vc_s; 
	};
	return F_XC_LDA<2>(inp, v);
}

double U_VWN5_c_E(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	// zeta = 0 constants
	const double A_0  = (1 - log(2)) / (M_PI * M_PI);
	const double x0_0 = -0.10498;
	const double b_0  =  3.72744;
	const double c_0  =  12.9352;
	const double X0_0 = x0_0 * x0_0 + b_0 * x0_0 + c_0;
	const double Q_0  = sqrt(4 * c_0 - b_0 * b_0);

	// zeta = 1 constants
	const double A_1  =  A_0 / 2;
	const double x0_1 = -0.32500;
	const double b_1  =  7.06042;
	const double c_1  =  18.0578;
	const double X0_1 = x0_1 * x0_1 + b_1 * x0_1 + c_1;
	const double Q_1  = sqrt(4 * c_1 - b_1 * b_1);
	
	auto e = [A_0, x0_0, b_0, c_0, X0_0, Q_0, A_1, x0_1, b_1, c_1, X0_1, Q_1](double rho_a, double rho_b) 
	{
		double rho = rho_a + rho_b;
		if(rho < 1e-20) {return 0.0;}
		double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));

		double zeta = ( rho_a - rho_b ) / rho;
		double zeta4 = zeta * zeta * zeta * zeta;
		double f = f_zeta(zeta);
		double ddf0 = 4.0 / ( 9.0 * ( cbrt(2) - 1 ) );
		double alpha = VWN_alpha(x);

		// evaluate energy densities for zeta = 0, 1
		double X_0  = x * x + b_0 * x + c_0;
		double X_1  = x * x + b_1 * x + c_1;

		double ec_0 = rho * A_0 * (
			log(x * x / X_0) + (2 * b_0 / Q_0) * (1 - (2 * x0_0 + b_0) * x0_0 / X0_0) * 
			atan(Q_0 / (2 * x + b_0)) - (b_0 * x0_0 / X0_0) * log((x - x0_0) * (x - x0_0) / X_0)
		);
		double ec_1 = rho * A_1 * (
			log(x * x / X_1) + (2 * b_1 / Q_1) * (1 - (2 * x0_1 + b_1) * x0_1 / X0_1) * 
			atan(Q_1 / (2 * x + b_1)) - (b_1 * x0_1 / X0_1) * log((x - x0_1) * (x - x0_1) / X_1)
		);

		return ec_0 + rho * alpha * (f / ddf0) * (1 - zeta4) + (ec_1 - ec_0) * f * zeta4;
	};
	return E_XC_LDA<1>(inp, e);	
}

XC_ret R_VWN5(const XC_inp& inp){
	Matrix null;
	return {R_VWN5_c(inp).F_XC_1 + R_Slater_X(inp).F_XC_1, null};
}

double R_VWN5_E(const XC_inp& inp){
	return R_VWN5_c_E(inp) + R_Slater_X_E(inp);
}

XC_ret U_VWN5(const XC_inp& inp){
	XC_ret fx = U_Slater_X(inp);
	XC_ret fc = U_VWN5_c(inp);
	return {fx.F_XC_1 + fc.F_XC_1, fx.F_XC_2 + fc.F_XC_2};
}

double U_VWN5_E(const XC_inp& inp){
	return U_VWN5_c_E(inp) + U_Slater_X_E(inp);
}

XC_ret R_PW92_c(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	// const double A  = 0.031091;
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double a1 = 0.21370;
	const double b1 = 7.5957;
	const double b2 = 3.5876;
	const double b3 = 1.6382;
	const double b4 = 0.49294;

	auto v = [A, a1, b1, b2, b3, b4](double rho){
		if(rho < 1e-20) {return 0.0;}
		double rs = cbrt(3 / (4 * M_PI * rho));
		double Q0  = -2 * A * (1 + a1 * rs);
		double Q1  =  2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs);
		double Q1p =      A * (b1 / sqrt(rs) + 2 * b2 + 3 * b3 * sqrt(rs) + 4 * b4 * rs);
		
		return - 2 * A * (1 + a1 * rs) * log(1 + 1 / (2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs))) 
			   - (rs / 3) * (-2 * A * a1 * log(1 + 1 / Q1) - Q0 * Q1p / (Q1 * Q1 + Q1));
	};
	return F_XC_LDA<0>(inp, v);
}

double R_PW92_c_E(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	//const double A  = 0.031091;
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double a1 = 0.21370;
	const double b1 = 7.5957;
	const double b2 = 3.5876;
	const double b3 = 1.6382;
	const double b4 = 0.49294;

	auto e = [A, a1, b1, b2, b3, b4](double rho){ 
		if(rho < 1e-20) {return 0.0;}
		double rs = cbrt(3 / (4 * M_PI * rho));
		return -rho * 2 * A * (1 + a1 * rs) * log(1 + 1 / (2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs)));
	};
	return E_XC_LDA<0>(inp, e);
}

XC_ret U_PW92_c(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	// zeta = 0
	// const double A_0  = 0.031091;
	const double A_0  = (1 - log(2)) / (M_PI * M_PI);
	const double a1_0 = 0.21370;
	const double b1_0 = 7.5957;
	const double b2_0 = 3.5876;
	const double b3_0 = 1.6382;
	const double b4_0 = 0.49294;
	// zeta = 1
	// const double A_1  = 0.015545;
	const double A_1  = A_0 / 2;
	const double a1_1 = 0.20548;
	const double b1_1 = 14.1189;
	const double b2_1 = 6.1977;
	const double b3_1 = 3.3662;
	const double b4_1 = 0.62517;

	auto v = [A_0, a1_0, b1_0, b2_0, b3_0, b4_0, A_1, a1_1, b1_1, b2_1, b3_1, b4_1] (double rho_a, double rho_b, int spin){
		double rho = rho_a + rho_b;
		if (rho < 1e-20) {return 0.0;}
		double rs = cbrt(3 / (4 * M_PI * rho));

		int sgn_spin;
		if(spin == 0)		{sgn_spin =  1;}
		else if (spin == 1) {sgn_spin = -1;}
		else {assert((spin == 0) || (spin == 1));}
		double zeta = (rho_a - rho_b) / rho;
		double zeta3 = zeta * zeta * zeta;
		double zeta4 = zeta3 * zeta;
		double f = f_zeta(zeta);
		double df = df_zeta(zeta);
		double ddf0 = 4.0 / ( 9.0 * ( cbrt(2) - 1 ) );
		double alpha = PW92_alpha(rs);
		double dalpha_drs = PW92_dalpha_drs(rs);

		double eps_0 = -2 * A_0 * (1 + a1_0 * rs) * log(1 + 1 / (2 * A_0 * 
					  (b1_0 * sqrt(rs) + b2_0 * rs + b3_0 * sqrt(intpow(rs, 3)) + b4_0 * rs * rs)));
		double eps_1 = -2 * A_1 * (1 + a1_1 * rs) * log(1 + 1 / (2 * A_1 * 
					  (b1_1 * sqrt(rs) + b2_1 * rs + b3_1 * sqrt(intpow(rs, 3)) + b4_1 * rs * rs)));
		double eps = eps_0 + alpha * (f / ddf0) * (1 - zeta4) + (eps_1 - eps_0) * f * zeta4;

		double Q0_0   = -2 * A_0 * (1 + a1_0 * rs);
		double Q1_0   =  2 * A_0 * (b1_0 * sqrt(rs) + b2_0 * rs + b3_0 * sqrt(intpow(rs, 3)) + b4_0 * rs * rs);
		double Q1p_0  =      A_0 * (b1_0 / sqrt(rs) + 2 * b2_0 + 3 * b3_0 * sqrt(rs) + 4 * b4_0 * rs);
		double deps_0 = -2 * A_0 * a1_0 * log(1 + 1 / Q1_0) - Q0_0 * Q1p_0 / (Q1_0 * Q1_0 + Q1_0);
		
		double Q0_1   = -2 * A_1 * (1 + a1_1 * rs);
		double Q1_1   =  2 * A_1 * (b1_1 * sqrt(rs) + b2_1 * rs + b3_1 * sqrt(intpow(rs, 3)) + b4_1 * rs * rs);
		double Q1p_1  =      A_1 * (b1_1 / sqrt(rs) + 2 * b2_1 + 3 * b3_1 * sqrt(rs) + 4 * b4_1 * rs);
		double deps_1 = -2 * A_1 * a1_1 * log(1 + 1 / Q1_1) - Q0_1 * Q1p_1 / (Q1_1 * Q1_1 + Q1_1);

		double deps_dr = deps_0 * (1 - f * zeta4) + deps_1 * f * zeta4 + dalpha_drs * (f / ddf0) * (1 - zeta4);
		double deps_dz = 4 * zeta3 * f * (eps_1 - eps_0 - alpha / ddf0) + df * (zeta4 * (eps_1 - eps_0) + (1 - zeta4) * alpha / ddf0);

		return eps - (rs / 3) * deps_dr - (zeta - sgn_spin) * deps_dz;
	};
	return F_XC_LDA<2>(inp, v);
}

double U_PW92_c_E(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	// zeta = 0
	// const double A_0  = 0.031091;
	const double A_0  = (1 - log(2)) / (M_PI * M_PI);
	const double a1_0 = 0.21370;
	const double b1_0 = 7.5957;
	const double b2_0 = 3.5876;
	const double b3_0 = 1.6382;
	const double b4_0 = 0.49294;
	// zeta = 1
	// const double A_1  = 0.015545;
	const double A_1  = A_0 / 2;
	const double a1_1 = 0.20548;
	const double b1_1 = 14.1189;
	const double b2_1 = 6.1977;
	const double b3_1 = 3.3662;
	const double b4_1 = 0.62517;

	auto e = [A_0, a1_0, b1_0, b2_0, b3_0, b4_0, A_1, a1_1, b1_1, b2_1, b3_1, b4_1] (double rho_a, double rho_b){
		double rho = rho_a + rho_b;
		if (rho < 1e-20) {return 0.0;}
		double rs = cbrt(3 / (4 * M_PI * rho));
		
		double zeta = (rho_a - rho_b) / rho;
		double zeta4 = zeta * zeta * zeta * zeta;
		double f = f_zeta(zeta);
		double ddf0 = 4.0 / ( 9.0 * ( cbrt(2) - 1 ) );
		double alpha = PW92_alpha(rs);

		double ec_0 = -rho * 2 * A_0 * (1 + a1_0 * rs) * log(1 + 1 / (2 * A_0 * 
					  (b1_0 * sqrt(rs) + b2_0 * rs + b3_0 * sqrt(intpow(rs, 3)) + b4_0 * rs * rs)));
		double ec_1 = -rho * 2 * A_1 * (1 + a1_1 * rs) * log(1 + 1 / (2 * A_1 * 
					  (b1_1 * sqrt(rs) + b2_1 * rs + b3_1 * sqrt(intpow(rs, 3)) + b4_1 * rs * rs)));

		return ec_0 + rho * alpha * (f / ddf0) * (1 - zeta4) + (ec_1 - ec_0) * f * zeta4;
	};
	return E_XC_LDA<1>(inp, e);
}

XC_ret R_PW92(const XC_inp& inp){
	Matrix null;
	return {R_PW92_c(inp).F_XC_1 + R_Slater_X(inp).F_XC_1, null};
}

double R_PW92_E(const XC_inp& inp){
	return R_PW92_c_E(inp) + R_Slater_X_E(inp);
}

XC_ret U_PW92(const XC_inp& inp){
	XC_ret fx = U_Slater_X(inp);
	XC_ret fc = U_PW92_c(inp);
	return {fx.F_XC_1 + fc.F_XC_1, fx.F_XC_2 + fc.F_XC_2};
}

double U_PW92_E(const XC_inp& inp){
	return U_PW92_c_E(inp) + U_Slater_X_E(inp);
}

// GGA ////////////////////////////////////////////////////////
XC_ret R_PBE_X(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	const double beta = 0.066725;
	const double kappa = 0.804;
	const double mu = beta * (M_PI * M_PI / 3);
	auto v = [kappa, mu](double rho, const std::vector<double>& grho, double phi1, double phi2, 
						 double gpx1, double gpy1, double gpz1,
					 	 double gpx2, double gpy2, double gpz2) 
	{
		if (rho < 1e-16) {return 0.0;}
		// Slater Exchange
		double v_LDA, e_LDA;
		v_LDA = -cbrt(3 * rho / M_PI);
		e_LDA = -(3.0 / 4.0) * cbrt(3.0 / M_PI) * cbrt(rho * rho * rho * rho);

		// Enhancement Factor
		double grho2, kF, s2;
		grho2 = grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2];
		kF = cbrt(3 * M_PI * M_PI * rho);
		s2 = grho2 / (4 * kF * kF * rho * rho);
		if (s2 < 1e-16) {return phi1 * v_LDA * phi2;}	

		double ds2_drho, ds2_dgrho2, FX_d, FX, dFX_ds2, dFX_drho, dFX_dgrho2;
		ds2_drho = -8.0 * s2 / (3.0 * rho);
		ds2_dgrho2 = s2 / grho2;	
		FX_d = 1 + mu * s2 / kappa;	
		FX = 1 + kappa - kappa / FX_d;
		dFX_ds2 = mu / (FX_d * FX_d);
		dFX_drho = dFX_ds2 * ds2_drho; 
		dFX_dgrho2 = dFX_ds2 * ds2_dgrho2;
 
		double de_drho = v_LDA * FX + e_LDA * dFX_drho;
		double de_dgrho2 = e_LDA * dFX_dgrho2;

		return (phi1 * de_drho * phi2) + 2 * de_dgrho2 * 
			   (phi1 * (grho[0] * gpx2 + grho[1] * gpy2 + grho[2] * gpz2) + 
			    phi2 * (grho[0] * gpx1 + grho[1] * gpy1 + grho[2] * gpz1));
	};
	return F_XC_GGA<0>(inp, v);
}

double R_PBE_X_E(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	const double beta = 0.066725;
	const double kappa = 0.804;
	const double mu = beta * (M_PI * M_PI / 3);
	auto e = [kappa, mu](double rho, std::vector<double> grho) 
	{
		if (rho < 1e-20) {return 0.0;}
		// Slater Exchange
		double e_LDA;
		e_LDA = -(3.0 / 4.0) * cbrt(3.0 / M_PI) * cbrt(rho * rho * rho * rho);

		// Enhancement Factor
		double grho2, kF, s2;
		grho2 = grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2];
		kF = cbrt(3 * M_PI * M_PI * rho);
		s2 = grho2 / (4 * kF * kF * rho * rho);

		return e_LDA * (1 + kappa - kappa / (1 + mu * s2 / kappa));
	};
	return E_XC_GGA<0>(inp, e);
}

XC_ret U_PBE_X(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	const double beta = 0.066725;
	const double kappa = 0.804;
	const double mu = beta * (M_PI * M_PI / 3);
	auto v = [kappa, mu](double rho_a, double rho_b, const std::vector<double>& grho_a, const std::vector<double>& grho_b, 
						 double phi1, double phi2, 
						 double gpx1, double gpy1, double gpz1,
					 	 double gpx2, double gpy2, double gpz2, int spin) 
	{
		double rho_s;
		std::vector<double> grho_s;
		if     (spin == 0) {rho_s = 2 * rho_a; grho_s = grho_a;}
		else if(spin == 1) {rho_s = 2 * rho_b; grho_s = grho_b;}
		else{assert((spin==0) || (spin==1));}
		if (rho_s / 2 < 1e-20) {return 0.0;}
		// Slater Exchange
		double v_LDA, e_LDA;
		v_LDA = -cbrt(3 * rho_s / M_PI);
		e_LDA = -(3.0 / 4.0) * cbrt(3.0 / M_PI) * cbrt(rho_s * rho_s * rho_s * rho_s);

		// Enhancement Factor
		double grho2, kF, s2;
		grho2 = 4 * (grho_s[0] * grho_s[0] + grho_s[1] * grho_s[1] + grho_s[2] * grho_s[2]);
		kF = cbrt(3 * M_PI * M_PI * rho_s);
		s2 = grho2 / (4 * kF * kF * rho_s * rho_s);

		double ds2_drho, ds2_dgrho2, FX_d, FX, dFX_ds2, dFX_drho, dFX_dgrho2;
		ds2_drho = -8.0 * s2 / (3.0 * rho_s);
		ds2_dgrho2 = s2 / grho2;	
		FX_d = 1 + mu * s2 / kappa;	
		FX = 1 + kappa - kappa / FX_d;
		dFX_ds2 = mu / (FX_d * FX_d);
		dFX_drho = dFX_ds2 * ds2_drho;
		dFX_dgrho2 = dFX_ds2 * ds2_dgrho2;
 
		double de_drho = 2 * (v_LDA * FX + e_LDA * dFX_drho);
		double de_dgrho2 = 4 * e_LDA * dFX_dgrho2;

		return 0.5 * (phi1 * de_drho * phi2 + 2 * de_dgrho2 * 
			   (phi1 * (grho_s[0] * gpx2 + grho_s[1] * gpy2 + grho_s[2] * gpz2) + 
			    phi2 * (grho_s[0] * gpx1 + grho_s[1] * gpy1 + grho_s[2] * gpz1)));
	};
	return F_XC_GGA<1>(inp, v);
}

double U_PBE_X_E(const XC_inp& inp){
	const double beta = 0.066725;
	const double kappa = 0.804;
	const double mu = beta * (M_PI * M_PI / 3);
	auto e = [kappa, mu](double rho_a, double rho_b, const std::vector<double>& grho_a, const std::vector<double>& grho_b) 
	{
		double rho = rho_a + rho_b;
		if (rho < 1e-20) {return 0.0;}
		// Slater Exchange
		double e_LDA_a = -(3.0 / 4.0) * cbrt(6.0 / M_PI) * cbrt(rho_a * rho_a * rho_a * rho_a);
		double e_LDA_b = -(3.0 / 4.0) * cbrt(6.0 / M_PI) * cbrt(rho_b * rho_b * rho_b * rho_b);

		// Enhancement Factor
		double grho2_a, grho2_b, kF_a, kF_b, s2_a, s2_b, FX_a, FX_b;
		grho2_a = 4 * (grho_a[0] * grho_a[0] + grho_a[1] * grho_a[1] + grho_a[2] * grho_a[2]);
		grho2_b = 4 * (grho_b[0] * grho_b[0] + grho_b[1] * grho_b[1] + grho_b[2] * grho_b[2]);
		kF_a = cbrt(3 * M_PI * M_PI * rho_a * 2);
		kF_b = cbrt(3 * M_PI * M_PI * rho_b * 2);
		s2_a = grho2_a / (4 * kF_a * kF_a * rho_a * 2 * rho_a * 2);
		s2_b = grho2_b / (4 * kF_b * kF_b * rho_b * 2 * rho_b * 2);

		FX_a = (1 + kappa - kappa / (1 + mu * s2_a / kappa));
		FX_b = (1 + kappa - kappa / (1 + mu * s2_b / kappa));
		return e_LDA_a * FX_a + e_LDA_b * FX_b;
	};
	return E_XC_GGA<1>(inp, e);
}

XC_ret R_PBE_c(const XC_inp& inp){
	// PW92 params
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double a1 = 0.21370;
	const double b1 = 7.5957;
	const double b2 = 3.5876;
	const double b3 = 1.6382;
	const double b4 = 0.49294;
	// PBE params
	const double beta = 0.066725;
	const double gamma = A;

	auto v = [A, a1, b1, b2, b3, b4, beta, gamma](double rho, const std::vector<double>& grho, double phi1, double phi2,
						 						  double gpx1, double gpy1, double gpz1,
					 	 						  double gpx2, double gpy2, double gpz2) 
	{
		if(rho < 1e-20) {return 0.0;}
		double rs, e_LDA, Q0, Q1, Q1p, v_LDA, deps_dn, kF, ks, grho2, t2, A_PBE, dnm, Q, H;
		rs = cbrt(3 / (4 * M_PI * rho));
		e_LDA = -rho * 2 * A * (1 + a1 * rs) * log(1 + 1 / (2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs)));

		Q0  = -2 * A * (1 + a1 * rs);
		Q1  =  2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs);
		Q1p =      A * (b1 / sqrt(rs) + 2 * b2 + 3 * b3 * sqrt(rs) + 4 * b4 * rs);
		
		v_LDA = - 2 * A * (1 + a1 * rs) * log(1 + 1 / (2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs))); 
		deps_dn = -(rs / 3) * (-2 * A * a1 * log(1 + 1 / Q1) - Q0 * Q1p / (Q1 * Q1 + Q1));
		v_LDA += deps_dn;

		kF = cbrt(3 * M_PI * M_PI * rho);
		ks = sqrt(4 * kF / M_PI);
		grho2 = grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2];
		t2 = grho2 / (4 * ks * ks * rho * rho);
		A_PBE = (beta / gamma) / (exp(-e_LDA / (rho * gamma)) - 1);
		if(A_PBE > 1e50) {return phi1 * v_LDA * phi2;}
		dnm = 1 + A_PBE * t2 + A_PBE * A_PBE * t2 * t2;
		Q = 1 + (beta / gamma) * t2 * (1 + A_PBE * t2) / dnm;
		H = gamma * log(Q);

		double dH_dt2, dt2_dn, dH_dA_PBE, dA_PBE_deps, dt2_dgrho2;
		dH_dt2 = (beta / Q) * (1 + 2 * A * t2) / (dnm * dnm);
		dt2_dn = -(7.0 / 3.0) * t2 / rho;
		dH_dA_PBE = -(beta / Q) * A_PBE * t2 * t2 * t2 * (2 + A_PBE * t2) / (dnm * dnm);
		dA_PBE_deps = A_PBE * (A_PBE + beta / gamma) / beta;
		dt2_dgrho2 = t2 / grho2;

		return phi1 * (v_LDA + H + rho * (dH_dt2 * dt2_dn + dH_dA_PBE * dA_PBE_deps * deps_dn)) * phi2 + 
			   2 * rho * dH_dt2 * dt2_dgrho2 * (phi2 * (grho[0] * gpx1 + grho[1] * gpy1 + grho[3] * gpz1) + 
			   phi1 * (grho[0] * gpx2 + grho[1] * gpy2 + grho[2] * gpz2));
	};
	return F_XC_GGA<0>(inp, v);
}

double R_PBE_c_E(const XC_inp& inp){
	// PW92 params
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double a1 = 0.21370;
	const double b1 = 7.5957;
	const double b2 = 3.5876;
	const double b3 = 1.6382;
	const double b4 = 0.49294;
	// PBE params
	const double beta = 0.066725;
	const double gamma = A;

	auto e = [A, a1, b1, b2, b3, b4, beta, gamma](double rho, const std::vector<double>& grho){
		if(rho < 1e-20) {return 0.0;}
		double rs, e_LDA, kF, ks, grho2, t2, A_PBE, Q;
		rs = cbrt(3 / (4 * M_PI * rho));
		e_LDA = -rho * 2 * A * (1 + a1 * rs) * log(1 + 1 / (2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs)));
		kF = cbrt(3 * M_PI * M_PI * rho);
		ks = sqrt(4 * kF / M_PI);
		grho2 = grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2];
		t2 = grho2 / (4 * ks * ks * rho * rho);
		A_PBE = (beta / gamma) / (exp(-e_LDA / (rho * gamma)) - 1);
		if(A_PBE > 1e50) {return e_LDA;}
		Q = 1 + (beta / gamma) * t2 * (1 + A_PBE * t2) / (1 + A_PBE * t2 + A_PBE * A_PBE * t2 * t2);
		return e_LDA + rho * gamma * log (Q);
	};
	return E_XC_GGA<0>(inp, e);
}

XC_ret R_PBE(const XC_inp& inp){
	Matrix null;
	return {R_PBE_c(inp).F_XC_1 + R_PBE_X(inp).F_XC_1, null};
}

double R_PBE_E(const XC_inp& inp){
	return R_PBE_c_E(inp) + R_PBE_X_E(inp);
}

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
