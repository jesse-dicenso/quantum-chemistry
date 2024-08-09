#include "1e.hpp"

double Sp(std::vector<int> L1, std::vector<int> L2, double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2){
	std::vector<double> QAB({xyz1[0]-xyz2[0], xyz1[1]-xyz2[1], xyz1[2]-xyz2[2]});
	return	E(L1[0], L2[0], 0, exp1, exp2, QAB[0]) * 
		E(L1[1], L2[1], 0, exp1, exp2, QAB[1]) *
		E(L1[2], L2[2], 0, exp1, exp2, QAB[2]) *
		pow((M_PI/(exp1+exp2)),1.5);
}

double S(GF g1, GF g2){
        double sum = 0;
        for(int i = 0; i < g1.exps.size(); i++){
                for(int j = 0; j < g2.exps.size(); j++){
			sum += g1.N[i] * g2.N[j] * g1.d[i] * g2.d[j] * Sp(g1.shell, g2.shell, g1.exps[i], g2.exps[j], g1.xyz, g2.xyz);
                }
        }
        return sum;
}

/*
double T(PGF A, PGF B){
	return 0;
}

double T(CGF A, CGF B){
        double sum = 0;
        for(int i = 0; i < A.len; i++){
                for(int j = 0; j < B.len; j++){
                        sum += A.d[i] * B.d[j] * (T(A.pgf[i], B.pgf[j]) / (A.pgf[i].N * B.pgf[j].N));
                }
        }
        return (A.N * B.N * sum);
}

double V(PGF A, PGF B){

	return 0;
}

double V(CGF A, CGF B){

	return 0;
}
*/
