#include "dft_helper.hpp"

double density(double x, double y, double z, const Molecule& mol, const Matrix& P){
	double rho = 0;
	double eval_gf;
	for(int i = 0; i < P.rows; i++){
			eval_gf = mol.AOs[i].evaluate(x, y, z);
			for(int j = i+1; j < P.cols; j++){
					rho += 2 * P.matrix[i][j] * eval_gf * mol.AOs[j].evaluate(x, y, z);
			}
			rho += P.matrix[i][i] * eval_gf * eval_gf;
	}
	return rho;
}

std::vector<double> density_gradient(double x, double y, double z, const Molecule& mol, const Matrix& P){
	std::vector<double> grad_rho = {0.0, 0.0, 0.0};
	std::vector<double> eval_grad_gf;
	double eval_gf;
	for(int i = 0; i < P.rows; i++){
			eval_grad_gf = mol.AOs[i].evaluate_gradient(x, y, z);
			for(int j = 0; j < P.cols; j++){
					eval_gf = P.matrix[i][j] * mol.AOs[j].evaluate(x, y, z);
					grad_rho[0] += eval_grad_gf[0] * eval_gf;
					grad_rho[1] += eval_grad_gf[1] * eval_gf;
					grad_rho[2] += eval_grad_gf[2] * eval_gf;
			}
	}
	grad_rho[0] *= 2;
	grad_rho[1] *= 2;
	grad_rho[2] *= 2;
	return grad_rho;
}

double f_zeta(double zeta){
	return ( cbrt(intpow(1 + zeta, 4)) + cbrt(intpow(1 - zeta, 4)) - 2.0 ) / ( cbrt(intpow(2.0, 4)) - 2.0 );
}

double df_zeta(double zeta){
	return 2 * ( cbrt(1 + zeta) - cbrt(1 - zeta) ) / (3 * (cbrt(2) - 1));
}

// VWN spin stiffness
double VWN_alpha(double x){
	const double A  = -1 / (6 * M_PI * M_PI);
	const double x0 = -0.00475840;
	const double b  =  1.13107;
	const double c  =  13.0045;
	const double Q  = sqrt(4 * c - b * b);
	const double X0 = x0 * x0 + b * x0 + c;	

	double X  = x * x + b * x + c;
	return A * (
		log(x * x / X) + (2 * b / Q) * (1 - (2 * x0 + b) * x0 / X0) * atan(Q / (2 * x + b)) - (b * x0 / X0) * log((x - x0) * (x - x0) / X)
	);
}

// VWN spin stiffness, derivative w.r.t. density
double VWN_dalpha_drho(double x, double n){	
	const double A  = -1 / (6 * M_PI * M_PI);
	const double x0 = -0.00475840;
	const double b  =  1.13107;
	const double c  =  13.0045;
	const double Q  = sqrt(4 * c - b * b);
	const double X0 = x0 * x0 + b * x0 + c;

	double X  = x * x + b * x + c;
	return (A / (3 * X * n)) * (b / ((1 / x0) - (1 / x)) - c);
}
