#include "scfgen.hpp"

Matrix overlap(std::vector<GF> phis){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			M.matrix[i][j] = S(phis[i],phis[j]);
		}
	}
	return M;
}

Matrix kinetic(std::vector<GF> phis){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			M.matrix[i][j] = T(phis[i],phis[j]);
		}
	}
	return M;
}

Matrix nuclear(std::vector<GF> phis, std::vector<int> Zvals, std::vector<std::vector<double>> xyzN){
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

double nucrepl(std::vector<int> Z, std::vector<std::vector<double>> xyzN){
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
