#include "dft_helper.hpp"

double density(double x, double y, double z, const std::vector<double>& phis, const Matrix& P){
	double rho = 0;
	for(int i = 0; i < P.rows; i++){
			for(int j = i+1; j < P.cols; j++){
					rho += 2 * P.matrix[i][j] * phis[i] * phis[j];
			}
			rho += P.matrix[i][i] * phis[i] * phis[i];
	}
	return rho;
}

std::vector<double> density_gradient(double x, double y, double z, const std::vector<double>& phis, 
									 const std::vector<double>& g_phis, const Matrix& P)
{
	double temp;
	std::vector<double> grad_rho = {0.0, 0.0, 0.0};
	for(int i = 0; i < P.rows; i++){
			for(int j = 0; j < P.cols; j++){
					temp = P.matrix[i][j] * phis[j];
					grad_rho[0] += g_phis[0] * temp;
					grad_rho[1] += g_phis[1] * temp;
					grad_rho[2] += g_phis[2] * temp;
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
	return -(A / (3 * X * n)) * (c - b * x0 * x / (x - x0));
}

// PW92 spin stiffness
double PW92_alpha(double rs){
	// const double A  = 0.016887;
	const double A  = 1 / (6 * M_PI * M_PI);
	const double a1 = 0.11125;
	const double b1 = 10.357;
	const double b2 = 3.6231;
	const double b3 = 0.88026;
	const double b4 = 0.49671;

	return 2 * A * (1 + a1 * rs) * log(1 + 1 / (2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs)));
}

// PW92 spin stiffness, derivative w.r.t. rs
double PW92_dalpha_drs(double rs){
	// const double A  = 0.016887;
	const double A  = 1 / (6 * M_PI * M_PI);
	const double a1 = 0.11125;
	const double b1 = 10.357;
	const double b2 = 3.6231;
	const double b3 = 0.88026;
	const double b4 = 0.49671;

	double Q0  = -2 * A * (1 + a1 * rs);
	double Q1  =  2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs);
	double Q1p =      A * (b1 / sqrt(rs) + 2 * b2 + 3 * b3 * sqrt(rs) + 4 * b4 * rs);
	return 2 * A * a1 * log(1 + 1 / Q1) + Q0 * Q1p / (Q1 * Q1 + Q1);
}

// Old density function
double density2(double x, double y, double z, const Molecule& mol, const Matrix& P){
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
