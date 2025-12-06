#include "xc.hpp"

///// Exchange only /////
Matrix R_HF_X(const Matrix& P, const std::vector<std::vector<std::vector<std::vector<double>>>>& g){
	Matrix K(P.rows, P.cols);
	double sum;
	for(int mu = 0; mu < K.rows; mu++){
		for(int nu = 0; nu < K.cols; nu++){
			sum = 0;
			for(int ld = 0; ld < K.rows; ld++){
				for(int sg = 0; sg < K.cols; sg++){
					sum -= P.matrix[ld][sg] * g[mu][ld][sg][nu];
				}
			}
			K.matrix[mu][nu] = 0.5 * sum;
		}
	}
	return K;
}

Matrix U_HF_X_s(const Matrix& Ps, const std::vector<std::vector<std::vector<std::vector<double>>>>& g){
	Matrix K(Ps.rows, Ps.cols);
	double sum;
	for(int mu = 0; mu < K.rows; mu++){
		for(int nu = 0; nu < K.cols; nu++){
			sum = 0;
			for(int ld = 0; ld < K.rows; ld++){
				for(int sg = 0; sg < K.cols; sg++){
					sum -= Ps.matrix[ld][sg] * g[mu][ld][sg][nu];
				}
			}
			K.matrix[mu][nu] = sum;
		}
	}
	return K;
}

Matrix R_Slater_X(const Matrix& P, const grid& g){
	Matrix K(P.rows, P.cols);
	//
	return K;
}

///// Correlation only /////

///// Exchange and Correlation /////
/*
	VWN5
	rs   = cbrt(3 / (4 * M_PI * rho));
	zeta = (rho_up - rho_down) / rho;
*/
