#include "1e.hpp"

double S(PGF A, PGF B){
	double p = A.exp + B.exp;
	std::vector<double> QAB({A.xyz[0]-B.xyz[0], A.xyz[1]-B.xyz[1], A.xyz[2]-B.xyz[2]});
	return pow((M_PI/(A.exp+B.exp)),1.5) * E(A.shell[0], B.shell[0], 0, A.exp, B.exp, QAB[0])
							 * E(A.shell[1], B.shell[1], 0, A.exp, B.exp, QAB[1]) 
							 * E(A.shell[2], B.shell[2], 0, A.exp, B.exp, QAB[2]);	
}

double S(CGF A, CGF B){
        double sum = 0;
        for(int i = 0; i < A.len; i++){
                for(int j = 0; j < B.len; j++){
                        sum += A.d[i] * B.d[j] * S(A.pgf[i], B.pgf[j]);
                }
        }
        return sum;
}


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
