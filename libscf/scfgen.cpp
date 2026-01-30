#include "scfgen.hpp"

Matrix overlap(const std::vector<GF>& phis){
	const int size_p = phis.size();
	Matrix M(size_p, size_p);
	for(int i = 0; i < size_p; i++){
		for(int j = 0; j < size_p; j++){
			M.matrix[i][j] = S(phis[i],phis[j]);
		}
	}
	return M;
}

Matrix kinetic(const std::vector<GF>& phis){
	const int size_p = phis.size();
	Matrix M(size_p, size_p);
	for(int i = 0; i < size_p; i++){
		for(int j = 0; j < size_p; j++){
			M.matrix[i][j] = T(phis[i],phis[j]);
		}
	}
	return M;
}

Matrix nuclear(const std::vector<GF>& phis, const std::vector<int>& Zvals, const std::vector<std::vector<double>>& xyzN){
	const int size_p = phis.size();
	const int size_z = Zvals.size();
	Matrix M(size_p, size_p);
	for(int i = 0; i < size_p; i++){
		for(int j = 0; j < size_p; j++){
			for(int k = 0; k < size_z; k++){
				M.matrix[i][j] += -Zvals[k]*V(phis[i],phis[j],xyzN[k]);
			}
		}
	}
	return M;
}

Matrix coulomb(const Matrix& P, const std::vector<std::vector<std::vector<std::vector<double>>>>& g){
	Matrix J(P.rows, P.cols);
	double sum;
	for(int mu = 0; mu < J.rows; mu++){
		for(int nu = 0; nu < J.cols; nu++){
			sum = 0;
			for(int ld = 0; ld < J.rows; ld++){
				for(int sg = 0; sg < J.cols; sg++){
					sum += P.matrix[ld][sg] * g[mu][nu][sg][ld];
				}
			}
			J.matrix[mu][nu] = sum;
		}
	}
	return J;
}

Matrix fock(const Matrix& Hcore, const Matrix& J, const Matrix& K){
	return Hcore + J + K;
}

double nucrepl(const std::vector<int>& Z, const std::vector<std::vector<double>>& xyzN){
	const int size_z = Z.size();
	double sum = 0;
	double Rij;
	std::vector<double> Ri;
	std::vector<double> Rj;
	for(int i = 0; i < size_z; i++){
		for(int j = (i+1); j < size_z; j++){
			Ri = xyzN[i];
			Rj = xyzN[j];
			Rij = sqrt((Ri[0]-Rj[0])*(Ri[0]-Rj[0]) + 
				   (Ri[1]-Rj[1])*(Ri[1]-Rj[1]) + 
				   (Ri[2]-Rj[2])*(Ri[2]-Rj[2]));
			sum += Z[i]*Z[j] / Rij;
		}
	}
	return sum;
}

double E0(const XC& xc, const Matrix& Hcore, const Matrix& J){
	double sum = 0;
	for(int i = 0; i < xc.P->rows; i++){
		for(int j = 0; j < xc.P->cols; j++){
			sum += xc.P->matrix[j][i]*(Hcore.matrix[i][j]+0.5*J.matrix[i][j]);
		}
	}
	return sum + xc.E_XC;
}
