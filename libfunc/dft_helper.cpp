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
	return tau;
	/* B97M-V does not have the 1/2 prefactor! */
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
	const double X  = x * x + b * x + c;

	return -(A / (3 * X * n)) * (c - b * x0 * x / (x - x0));
}

// PW92 spin stiffness
double PW92_alpha(double rs){
	// const double A  = 0.016887;
	constexpr double A  = 1 / (6 * M_PI * M_PI);
	constexpr double a1 = 0.11125;
	constexpr double b1 = 10.357;
	constexpr double b2 = 3.6231;
	constexpr double b3 = 0.88026;
	constexpr double b4 = 0.49671;

	return 2 * A * (1 + a1 * rs) * log(1 + 1 / (2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs)));
}

// PW92 spin stiffness, derivative w.r.t. rs
double PW92_dalpha_drs(double rs){
	// const double A  = 0.016887;
	const double A  = 1 / (6 * M_PI * M_PI);
	constexpr double a1 = 0.11125;
	constexpr double b1 = 10.357;
	constexpr double b2 = 3.6231;
	constexpr double b3 = 0.88026;
	constexpr double b4 = 0.49671;

	const double Q0  = -2 * A * (1 + a1 * rs);
	const double Q1  =  2 * A * (b1 * sqrt(rs) + b2 * rs + b3 * sqrt(intpow(rs, 3)) + b4 * rs * rs);
	const double Q1p =      A * (b1 / sqrt(rs) + 2 * b2 + 3 * b3 * sqrt(rs) + 4 * b4 * rs);
	return 2 * A * a1 * log(1 + 1 / Q1) + Q0 * Q1p / (Q1 * Q1 + Q1);
}

// per particle!
double eps_c_pw92(double rho_a, double rho_b){
	// zeta = 0
	const double A_0  = (1 - log(2)) / (M_PI * M_PI);
	constexpr double a1_0 = 0.21370;
	constexpr double b1_0 = 7.5957;
	constexpr double b2_0 = 3.5876;
	constexpr double b3_0 = 1.6382;
	constexpr double b4_0 = 0.49294;
	// zeta = 1
	const double A_1  = A_0 / 2;
	constexpr double a1_1 = 0.20548;
	constexpr double b1_1 = 14.1189;
	constexpr double b2_1 = 6.1977;
	constexpr double b3_1 = 3.3662;
	constexpr double b4_1 = 0.62517;

	const double rho = rho_a + rho_b;
	if (rho < 1e-20) {return 0.0;}
	const double rs = cbrt(3 / (4 * M_PI * rho));
		
	const double zeta = (rho_a - rho_b) / rho;
	const double zeta4 = zeta * zeta * zeta * zeta;
	const double f = f_zeta(zeta);
	const double ddf0 = 4.0 / ( 9.0 * ( cbrt(2) - 1 ) );
	const double alpha = PW92_alpha(rs);

	const double epsc_0 = -2 * A_0 * (1 + a1_0 * rs) * log(1 + 1 / (2 * A_0 * 
				    (b1_0 * sqrt(rs) + b2_0 * rs + b3_0 * sqrt(intpow(rs, 3)) + b4_0 * rs * rs)));
	const double epsc_1 = -2 * A_1 * (1 + a1_1 * rs) * log(1 + 1 / (2 * A_1 * 
				  (b1_1 * sqrt(rs) + b2_1 * rs + b3_1 * sqrt(intpow(rs, 3)) + b4_1 * rs * rs)));

	return epsc_0 + alpha * (f / ddf0) * (1 - zeta4) + (epsc_1 - epsc_0) * f * zeta4;
}

double deps_c_dns_pw92(double rho_a, double rho_b, int spin){
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

	const double rho = rho_a + rho_b;
	if (rho < 1e-20) {return 0.0;}
	const double rs = cbrt(3 / (4 * M_PI * rho));

	// spin == 0 : take derivative w.r.t. rho_a
	// spin == 1 : take derivative w.r.t. rho_b
	const int sgn_spin = (spin==0 ? 1 : -1);
	const double zeta = (rho_a - rho_b) / rho;
	const double zeta3 = zeta * zeta * zeta;
	const double zeta4 = zeta3 * zeta;
	const double f = f_zeta(zeta);
	const double df = df_zeta(zeta);
	const double ddf0 = 4.0 / ( 9.0 * ( cbrt(2) - 1 ) );
	const double alpha = PW92_alpha(rs);
	const double dalpha_drs = PW92_dalpha_drs(rs);

	const double eps_0 = -2 * A_0 * (1 + a1_0 * rs) * log(1 + 1 / (2 * A_0 * 
				  (b1_0 * sqrt(rs) + b2_0 * rs + b3_0 * sqrt(intpow(rs, 3)) + b4_0 * rs * rs)));
	const double eps_1 = -2 * A_1 * (1 + a1_1 * rs) * log(1 + 1 / (2 * A_1 * 
				  (b1_1 * sqrt(rs) + b2_1 * rs + b3_1 * sqrt(intpow(rs, 3)) + b4_1 * rs * rs)));

	const double Q0_0   = -2 * A_0 * (1 + a1_0 * rs);
	const double Q1_0   =  2 * A_0 * (b1_0 * sqrt(rs) + b2_0 * rs + b3_0 * sqrt(intpow(rs, 3)) + b4_0 * rs * rs);
	const double Q1p_0  =      A_0 * (b1_0 / sqrt(rs) + 2 * b2_0 + 3 * b3_0 * sqrt(rs) + 4 * b4_0 * rs);
	const double deps_0 = -2 * A_0 * a1_0 * log(1 + 1 / Q1_0) - Q0_0 * Q1p_0 / (Q1_0 * Q1_0 + Q1_0);
		
	const double Q0_1   = -2 * A_1 * (1 + a1_1 * rs);
	const double Q1_1   =  2 * A_1 * (b1_1 * sqrt(rs) + b2_1 * rs + b3_1 * sqrt(intpow(rs, 3)) + b4_1 * rs * rs);
	const double Q1p_1  =      A_1 * (b1_1 / sqrt(rs) + 2 * b2_1 + 3 * b3_1 * sqrt(rs) + 4 * b4_1 * rs);
	const double deps_1 = -2 * A_1 * a1_1 * log(1 + 1 / Q1_1) - Q0_1 * Q1p_1 / (Q1_1 * Q1_1 + Q1_1);

	const double deps_dr = deps_0 * (1 - f * zeta4) + deps_1 * f * zeta4 + dalpha_drs * (f / ddf0) * (1 - zeta4);
	const double deps_dz = 4 * zeta3 * f * (eps_1 - eps_0 - alpha / ddf0) + df * (zeta4 * (eps_1 - eps_0) + (1 - zeta4) * alpha / ddf0);

	return -((rs / 3) * deps_dr + (zeta - sgn_spin) * deps_dz) / rho;
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
