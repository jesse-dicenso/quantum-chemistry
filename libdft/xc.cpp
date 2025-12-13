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
	{ "U_Slater", U_Slater_X }/*,
	{ "R_VWN5_c", R_VWN5_c },
	{ "U_VWN5_c", U_VWN5_c },
	{ "R_VWN5, R_VWN5" },
	{ "U_VWN5, U_VWN5" },
	*/
};

XC_ret F_XC(XC_inp* inp){
	assert(inp!=nullptr);
	return xc_v_register[inp->method](*inp);
}

std::unordered_map<std::string, std::function<double(const XC_inp&)>> xc_E_register = 
{
	// LDA //
	{ "R_Slater", R_Slater_X_E },
	{ "U_Slater", U_Slater_X_E }/*,
	{ "R_VWN5_c", R_VWN5_c_E },
	{ "U_VWN5_c", U_VWN5_c_E },
	{ "R_VWN5, R_VWN5_E" },
	{ "U_VWN5, U_VWN5_E" },
	*/
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
/*
	VWN5
	rs   = cbrt(3 / (4 * M_PI * rho));
	zeta = (rho_up - rho_down) / rho;
*/
