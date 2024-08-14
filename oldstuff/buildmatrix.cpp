#include "buildmatrix.hpp"

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
