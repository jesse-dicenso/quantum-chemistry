#include "xc.hpp"

XC_inp::XC_inp(const std::string& method_name){
	method = method_name;

	is_HF = is_LDA = is_GGA = false;
	assert((method.substr(0,2)=="R_") || (method.substr(0,2)=="U_"));
	if      (method.substr(2)=="HF"      )  {is_HF  = true;}
	else if((method.substr(2)=="Slater") || 
			(method.substr(2)=="VWN5_c"  ) || 
			(method.substr(2)=="VWN5"    )) {is_LDA = true;}
	else{
		std::cerr << "Error: method " << method << " not found!" << std::endl;
		assert(false);
	}
}

std::unordered_map<std::string, std::function<XC_ret(const XC_inp&)>> xc_v_register = 
{
	// HF //
	{ "R_HF", R_HF_X },
	{ "U_HF", U_HF_X },
	// LDA //
	{ "R_Slater", R_Slater_X },
	{ "U_Slater", U_Slater_X },
	{ "R_VWN5_c", R_VWN5_c },
	{ "U_VWN5_c", U_VWN5_c },
	{ "R_VWN5", R_VWN5 },
	{ "U_VWN5", U_VWN5 },
};

XC_ret F_XC(XC_inp* inp){
	assert(inp!=nullptr);
	return xc_v_register[inp->method](*inp);
}

std::unordered_map<std::string, std::function<double(const XC_inp&)>> xc_E_register = 
{
	// LDA //
	{ "R_Slater", R_Slater_X_E },
	{ "U_Slater", U_Slater_X_E },
	{ "R_VWN5_c", R_VWN5_c_E },
	{ "U_VWN5_c", U_VWN5_c_E },
	{ "R_VWN5", R_VWN5_E },
	{ "U_VWN5", U_VWN5_E },
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
	Matrix F_XC(inp.PT->rows, inp.PT->cols);
	Matrix null;
	auto integrand = [](double x, double y, double z, const Molecule& m, const Matrix& p, int idx1, int idx2) {
		return -m.AOs[idx1].evaluate(x,y,z) * cbrt(3 * density(x, y, z, m, p) / M_PI) * m.AOs[idx2].evaluate(x,y,z);
	};
	for(int i = 0; i < F_XC.rows; i++){
		F_XC.matrix[i][i] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PT, i, i);
		for(int j = 0; j < i; j++){
			F_XC.matrix[i][j] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PT, i, j);
			F_XC.matrix[j][i] = F_XC.matrix[i][j];
		}
	}
	return {F_XC, null};
}

double R_Slater_X_E(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	auto integrand = [](double x, double y, double z, const Molecule& m, const Matrix& p) {
		double rho = density(x, y, z, m, p);
		return cbrt(rho * rho * rho * rho);
	};
	return -(3.0/4.0) * cbrt(3.0 / M_PI) * integrate_quad(*inp.g, integrand, *inp.mol, *inp.PT);
}

XC_ret U_Slater_X(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	Matrix F_XC_a(inp.PA->rows, inp.PA->cols);
	Matrix F_XC_b(inp.PA->rows, inp.PB->cols);
	auto integrand = [](double x, double y, double z, const Molecule& m, const Matrix& p, int idx1, int idx2) {
		return -m.AOs[idx1].evaluate(x,y,z) * cbrt(6 * density(x, y, z, m, p) / M_PI) * m.AOs[idx2].evaluate(x,y,z);
	};
	for(int i = 0; i < F_XC_a.rows; i++){
		F_XC_a.matrix[i][i] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, i, i);
		F_XC_b.matrix[i][i] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PB, i, i);
		for(int j = 0; j < i; j++){
			F_XC_a.matrix[i][j] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, i, j);
			F_XC_b.matrix[i][j] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PB, i, j);
			F_XC_a.matrix[j][i] = F_XC_a.matrix[i][j];
			F_XC_b.matrix[j][i] = F_XC_b.matrix[i][j];
		}
	}
	return {F_XC_a, F_XC_b};
}

double U_Slater_X_E(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	auto integrand = [](double x, double y, double z, const Molecule& m, const Matrix& pa, const Matrix& pb) {
		double rho_a = density(x, y, z, m, pa);
		double rho_b = density(x, y, z, m, pb);
		return cbrt(rho_a * rho_a * rho_a * rho_a) + cbrt(rho_b * rho_b * rho_b * rho_b);
	};
	return -(3.0/4.0) * cbrt(6.0 / M_PI) * integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, *inp.PB);
}

XC_ret R_VWN5_c(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	Matrix F_XC(inp.PT->rows, inp.PT->cols);
	Matrix null;
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double x0 = -0.10498;
	const double b  =  3.72744;
	const double c  =  12.9352;
	const double X0 = x0 * x0 + b * x0 + c;
	const double Q  = sqrt(4 * c - b * b);
	auto integrand = [A, x0, b, c, X0, Q](double r_x, double r_y, double r_z, const Molecule& m, const Matrix& p, int idx1, int idx2) {
		double rho = density(r_x, r_y, r_z, m, p);
		if(rho < 1e-16) {return 0.0;}
		double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		double X  = x * x + b * x + c;
		double vc = (
			log(x * x / X) + (2 * b / Q) * (1 - (2 * x0 + b) * x0 / X0) * atan(Q / (2 * x + b)) - (b * x0 / X0) * log((x - x0) * (x - x0) / X)
		);
		vc -= (x / (3 * X)) * (c / x - b * x0 / (x - x0));
		vc *= A;
		return m.AOs[idx1].evaluate(r_x, r_y, r_z) * vc * m.AOs[idx2].evaluate(r_x, r_y, r_z); 
	};
	for(int i = 0; i < F_XC.rows; i++){
		F_XC.matrix[i][i] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PT, i, i);
		for(int j = 0; j < i; j++){
			F_XC.matrix[i][j] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PT, i, j);
			F_XC.matrix[j][i] = F_XC.matrix[i][j];
		}
	}
	return {F_XC, null};
}

double R_VWN5_c_E(const XC_inp& inp){
	assert((inp.PT!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	const double A  = (1 - log(2)) / (M_PI * M_PI);
	const double x0 = -0.10498;
	const double b  =  3.72744;
	const double c  =  12.9352;
	const double Q  = sqrt(4 * c - b * b);
	const double X0 = x0 * x0 + b * x0 + c;
	auto integrand = [A, x0, b, c, X0, Q](double r_x, double r_y, double r_z, const Molecule& m, const Matrix& p) {
		double rho = density(r_x, r_y, r_z, m, p);
		if(rho < 1e-16) {return 0.0;}
		double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		double X  = x * x + b * x + c;
		double ec = rho * A * (
			log(x * x / X) + (2 * b / Q) * (1 - (2 * x0 + b) * x0 / X0) * atan(Q / (2 * x + b)) - (b * x0 / X0) * log((x - x0) * (x - x0) / X)
		);
		return ec; 
	};
	return integrate_quad(*inp.g, integrand, *inp.mol, *inp.PT);
}

XC_ret U_VWN5_c(const XC_inp& inp){
	assert((inp.PA!=nullptr) && (inp.PB!=nullptr) && (inp.mol!=nullptr) && (inp.g!=nullptr));
	Matrix F_XC_a(inp.PA->rows, inp.PA->cols);
	Matrix F_XC_b(inp.PA->rows, inp.PB->cols);

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
	
	auto integrand = [A_0, x0_0, b_0, c_0, X0_0, Q_0, A_1, x0_1, b_1, c_1, X0_1, Q_1](double r_x, double r_y, double r_z, const Molecule& m, 
																			const Matrix& pa, const Matrix& pb, int idx1, int idx2, int spin) 
	{
		double rho_a = density(r_x, r_y, r_z, m, pa);
		double rho_b = density(r_x, r_y, r_z, m, pb);
		double rho = rho_a + rho_b;
		if(rho < 1e-16) {return 0.0;}
		double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));

		double zeta = ( rho_a - rho_b ) / rho;
		double zeta3 = zeta * zeta * zeta;
		// spin = 0 -> alpha, spin = 1 -> beta
		double dzeta_drho;
		if(spin == 0) {dzeta_drho = 2 * rho_b / (rho * rho) ;}
		else if(spin == 1) {dzeta_drho = -2 * rho_a / (rho * rho) ;}
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
					  
		return m.AOs[idx1].evaluate(r_x, r_y, r_z) * vc_s * m.AOs[idx2].evaluate(r_x, r_y, r_z); 
	};

	for(int i = 0; i < F_XC_a.rows; i++){
		F_XC_a.matrix[i][i] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, *inp.PB, i, i, 0);
		F_XC_b.matrix[i][i] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, *inp.PB, i, i, 1);
		for(int j = 0; j < i; j++){
			F_XC_a.matrix[i][j] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, *inp.PB, i, j, 0);
			F_XC_b.matrix[i][j] = integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, *inp.PB, i, j, 1);
			F_XC_a.matrix[j][i] = F_XC_a.matrix[i][j];
			F_XC_b.matrix[j][i] = F_XC_b.matrix[i][j];
		}
	}
	return {F_XC_a, F_XC_b};
}

double U_VWN5_c_E(const XC_inp& inp){
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
	
	auto integrand = [A_0, x0_0, b_0, c_0, X0_0, Q_0, A_1, x0_1, b_1, c_1, X0_1, Q_1](double r_x, double r_y, double r_z, const Molecule& m, 
																					  const Matrix& pa, const Matrix& pb) 
	{
		double rho_a = density(r_x, r_y, r_z, m, pa);
		double rho_b = density(r_x, r_y, r_z, m, pb);
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
	return integrate_quad(*inp.g, integrand, *inp.mol, *inp.PA, *inp.PB);	
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

// GGA //
