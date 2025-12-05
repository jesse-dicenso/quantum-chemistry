#include "scfgen.hpp"

Matrix overlap(const std::vector<GF>& phis){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			M.matrix[i][j] = S(phis[i],phis[j]);
		}
	}
	return M;
}

Matrix kinetic(const std::vector<GF>& phis){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			M.matrix[i][j] = T(phis[i],phis[j]);
		}
	}
	return M;
}

Matrix nuclear(const std::vector<GF>& phis, const std::vector<int>& Zvals, const std::vector<std::vector<double>>& xyzN){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			for(int k = 0; k < Zvals.size(); k++){
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
	double sum = 0;
	double Rij;
	std::vector<double> Ri;
	std::vector<double> Rj;
	for(int i = 0; i < Z.size(); i++){
		for(int j = (i+1); j < Z.size(); j++){
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
