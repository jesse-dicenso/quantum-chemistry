#include "dft_helper.hpp"

double density(const std::vector<double>& phis, const Matrix& P){
	double rho = 0;
	for(int i = 0; i < P.rows; i++){
			for(int j = i+1; j < P.cols; j++){
					rho += 2 * P.matrix[i][j] * phis[i] * phis[j];
			}
			rho += P.matrix[i][i] * phis[i] * phis[i];
	}
	return rho;
}

std::vector<double> density_gradient(const std::vector<double>& phis, const std::vector<double>& gpx, const std::vector<double>& gpy, 
									 const std::vector<double>& gpz, const Matrix& P)
{
	double temp;
	std::vector<double> grad_rho = {0.0, 0.0, 0.0};
	for(int i = 0; i < P.rows; i++){
			for(int j = 0; j < P.cols; j++){
					temp = P.matrix[i][j] * phis[j];
					grad_rho[0] += gpx[i] * temp;
					grad_rho[1] += gpy[i] * temp;
					grad_rho[2] += gpz[i] * temp;
			}
	}
	grad_rho[0] *= 2;
	grad_rho[1] *= 2;
	grad_rho[2] *= 2;
	return grad_rho;
}

double ke_density(const std::vector<double>& gpx, const std::vector<double>& gpy, const std::vector<double>& gpz, const Matrix& P){
	double tau = 0;
	for(int i = 0; i < P.rows; i++){
		for(int j = 0; j < P.cols; j++){
			tau += P.matrix[i][j] * (gpx[i] * gpx[j] + gpy[i] * gpy[j] + gpz[i] * gpz[j]);
		}
	}
	return 0.5 * tau;
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

double ke_density_ueg(double rho){
	return (3.0 / 5.0) * cbrt(36 * intpow(M_PI, 4) * intpow(rho, 5));
}

double e_X_ueg(double rho){
	return -(3.0 / 2.0) * cbrt(intpow(rho, 4) * 3.0 / (4.0 * M_PI));
}

// per particle!
double eps_c_pw92(double rho_a, double rho_b){
	// zeta = 0
	const double A_0  = (1 - log(2)) / (M_PI * M_PI);
	const double a1_0 = 0.21370;
	const double b1_0 = 7.5957;
	const double b2_0 = 3.5876;
	const double b3_0 = 1.6382;
	const double b4_0 = 0.49294;
	// zeta = 1
	const double A_1  = A_0 / 2;
	const double a1_1 = 0.20548;
	const double b1_1 = 14.1189;
	const double b2_1 = 6.1977;
	const double b3_1 = 3.3662;
	const double b4_1 = 0.62517;

	double rho = rho_a + rho_b;
	if (rho < 1e-20) {return 0.0;}
	double rs = cbrt(3 / (4 * M_PI * rho));
		
	double zeta = (rho_a - rho_b) / rho;
	double zeta4 = zeta * zeta * zeta * zeta;
	double f = f_zeta(zeta);
	double ddf0 = 4.0 / ( 9.0 * ( cbrt(2) - 1 ) );
	double alpha = PW92_alpha(rs);

	double epsc_0 = -2 * A_0 * (1 + a1_0 * rs) * log(1 + 1 / (2 * A_0 * 
				    (b1_0 * sqrt(rs) + b2_0 * rs + b3_0 * sqrt(intpow(rs, 3)) + b4_0 * rs * rs)));
	double epsc_1 = -2 * A_1 * (1 + a1_1 * rs) * log(1 + 1 / (2 * A_1 * 
				  (b1_1 * sqrt(rs) + b2_1 * rs + b3_1 * sqrt(intpow(rs, 3)) + b4_1 * rs * rs)));

	return epsc_0 + alpha * (f / ddf0) * (1 - zeta4) + (epsc_1 - epsc_0) * f * zeta4;
}

// Capital Phi in the VV10 paper
double VV10_kernel(double b, double C, double R2, double rho_1, double rho_2, double nrm_grho_1, double nrm_grho_2){
	double omega_p2_1 = 4 * M_PI * rho_1;
	double omega_p2_2 = 4 * M_PI * rho_2;
	double omega_g2_1 = C * intpow(nrm_grho_1 / rho_1, 4);
	double omega_g2_2 = C * intpow(nrm_grho_2 / rho_2, 4);
	double omega_0_1 = sqrt(omega_g2_1 + omega_p2_1 / 3);
	double omega_0_2 = sqrt(omega_g2_2 + omega_p2_2 / 3);
	double kappa_1 = b * intpow(cbrt(3 * M_PI * M_PI * rho_1), 2) / sqrt(omega_p2_1);
	double kappa_2 = b * intpow(cbrt(3 * M_PI * M_PI * rho_2), 2) / sqrt(omega_p2_2);

	double g_1 = omega_0_1 * R2 + kappa_1;
	double g_2 = omega_0_2 * R2 + kappa_2;

	return -3 / (2 * g_1 * g_2 * (g_1 + g_2));
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
