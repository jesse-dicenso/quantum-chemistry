#include "xc.hpp"

XC_inp::XC_inp(const std::string& method_name){
	method = method_name;

	is_HF = is_LDA = is_GGA = false;
	assert((method.substr(0,2)=="R_") || (method.substr(0,2)=="U_"));
	if      (method.substr(2)=="HF"     )  {is_HF  = true;}
	else if((method.substr(2)=="Slater" ) || 
			(method.substr(2)=="VWN5_c" ) || 
			(method.substr(2)=="VWN5"   ) || 
			(method.substr(2)=="PW92_c" ) || 
			(method.substr(2)=="PW92"   )) {is_LDA = true;}
	else{
		std::cerr << "Error: method " << method << " not found!" << std::endl;
		assert(false);
	}
}

std::unordered_map<std::string, std::function<XC_ret(const XC_inp&)>> xc_v_register = 
{
	{ "R_HF", R_HF_X },
	{ "U_HF", U_HF_X },
	{ "R_Slater", R_Slater_X },
	{ "U_Slater", U_Slater_X },
	{ "R_VWN5_c", R_VWN5_c },
	{ "U_VWN5_c", U_VWN5_c },
	{ "R_VWN5", R_VWN5 },
	{ "U_VWN5", U_VWN5 },
	{ "R_PW92_c", R_VWN5 },
	{ "U_PW92_c", R_VWN5 },
	{ "R_PW92", R_VWN5 },
	{ "U_PW92", R_VWN5 },
};

XC_ret F_XC(XC_inp* inp){
	assert(inp!=nullptr);
	return xc_v_register[inp->method](*inp);
}

std::unordered_map<std::string, std::function<double(const XC_inp&)>> xc_E_register = 
{
	{ "R_Slater", R_Slater_X_E },
	{ "U_Slater", U_Slater_X_E },
	{ "R_VWN5_c", R_VWN5_c_E },
	{ "U_VWN5_c", U_VWN5_c_E },
	{ "R_VWN5", R_VWN5_E },
	{ "U_VWN5", U_VWN5_E },
	{ "R_PW92_c", R_VWN5_E },
	{ "U_PW92_c", U_VWN5_E },
	{ "R_PW92", R_VWN5_E },
	{ "U_PW92", U_VWN5_E },
};

double E_XC(XC_inp* inp){
	assert(inp!=nullptr);
	return xc_E_register[inp->method](*inp);
}

// HF //
XC_ret R_HF_X(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.eris!=nullptr));
	Matrix F_XC(inp.PT->rows, inp.PT->cols);
	Matrix null;
	double sum;
	for(int mu = 0; mu < F_XC.rows; mu++){
		for(int nu = 0; nu < F_XC.cols; nu++){
			sum = 0;
			for(int ld = 0; ld < F_XC.rows; ld++){
				for(int sg = 0; sg < F_XC.cols; sg++){
					sum -= inp.PT->matrix[ld][sg] * (*inp.eris)[mu][ld][sg][nu];
				}
			}
			F_XC.matrix[mu][nu] = 0.5 * sum;
		}
	}
	return {F_XC, null};
}

XC_ret U_HF_X(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.eris!=nullptr));
	Matrix F_XC_a(inp.PA->rows, inp.PA->cols);
	Matrix F_XC_b(inp.PB->rows, inp.PB->cols);
	double sum_a, sum_b;
	for(int mu = 0; mu < F_XC_a.rows; mu++){
		for(int nu = 0; nu < F_XC_a.cols; nu++){
			sum_a = 0;
			sum_b = 0;
			for(int ld = 0; ld < F_XC_a.rows; ld++){
				for(int sg = 0; sg < F_XC_a.cols; sg++){
					sum_a -= inp.PA->matrix[ld][sg] * (*inp.eris)[mu][ld][sg][nu];
					sum_b -= inp.PB->matrix[ld][sg] * (*inp.eris)[mu][ld][sg][nu];
				}
			}
			F_XC_a.matrix[mu][nu] = sum_a;
			F_XC_b.matrix[mu][nu] = sum_b;
		}
	}
	return {F_XC_a, F_XC_b};
}

// LDA //
XC_ret R_Slater_X(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	auto v = [](double rho) {
		return -cbrt(3 * rho / M_PI);
	};
	return F_XC_LDA<0>(inp, v);
}

double R_Slater_X_E(const XC_inp& inp){
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
		if(rho < 1e-16) {return 0.0;}
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
		if(rho < 1e-16) {return 0.0;}
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
		if(rho < 1e-16) {return 0.0;}
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
		if(rho < 1e-16) {return 0.0;}
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
		if(rho < 1e-16) {return 0.0;}
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
		if(rho < 1e-16) {return 0.0;}
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
		if (rho < 1e-16) {return 0.0;}
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

	auto v = [A_0, a1_0, b1_0, b2_0, b3_0, b4_0, A_1, a1_1, b1_1, b2_1, b3_1, b4_1] (double rho_a, double rho_b){
		double rho = rho_a + rho_b;
		if (rho < 1e-16) {return 0.0;}
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
	return E_XC_LDA<1>(inp, v);
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
	return U_PW92_c_E(inp) + R_Slater_X_E(inp);
}

// GGA //
