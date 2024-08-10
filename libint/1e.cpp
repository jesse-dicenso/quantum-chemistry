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

double Tp(std::vector<int> L1, std::vector<int> L2, double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2){
	double Ix = exp2*(2*L2[0]+1)*Sp(L1,L2,exp1,exp2,xyz1,xyz2) 
		  - 2*exp2*exp2*Sp(L1,{(L2[0]+2),L2[1],L2[2]},exp1,exp2,xyz1,xyz2) 
		  - 0.5*L2[0]*(L2[0]-1)*Sp(L1,{(L2[0]-2+abs(L2[0]-2))/2,L2[1],L2[2]},exp1,exp2,xyz1,xyz2);
	double Iy = exp2*(2*L2[1]+1)*Sp(L1,L2,exp1,exp2,xyz1,xyz2) 
		  - 2*exp2*exp2*Sp(L1,{L2[0],(L2[1]+2),L2[2]},exp1,exp2,xyz1,xyz2) 
		  - 0.5*L2[1]*(L2[1]-1)*Sp(L1,{L2[0],(L2[1]-2+abs(L2[1]-2))/2,L2[2]},exp1,exp2,xyz1,xyz2);
	double Iz = exp2*(2*L2[2]+1)*Sp(L1,L2,exp1,exp2,xyz1,xyz2) 
		  - 2*exp2*exp2*Sp(L1,{L2[0],L2[1],(L2[2]+2)},exp1,exp2,xyz1,xyz2) 
		  - 0.5*L2[2]*(L2[2]-1)*Sp(L1,{L2[0],L2[1],(L2[2]-2+abs(L2[2]-2))/2},exp1,exp2,xyz1,xyz2);
	return (Ix + Iy + Iz);
}

double T(GF g1, GF g2){
        double sum = 0;
        for(int i = 0; i < g1.exps.size(); i++){
                for(int j = 0; j < g2.exps.size(); j++){
			sum += g1.N[i] * g2.N[j] * g1.d[i] * g2.d[j] * Tp(g1.shell, g2.shell, g1.exps[i], g2.exps[j], g1.xyz, g2.xyz);
                }
        }
        return sum;
}
/*
double V(GF g1, GF g2){

	return 0;
}

double V(GF g1, GF g2){

	return 0;
}
*/
