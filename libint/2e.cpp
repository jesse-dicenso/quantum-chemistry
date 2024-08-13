#include "2e.hpp"

double Gp(std::vector<int> L1, std::vector<int> L2, std::vector<int> L3, std::vector<int> L4, double exp1, double exp2, double exp3, double exp4, std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyz3, std::vector<double> xyz4){
	double p = exp1 + exp2;
	double q = exp3 + exp4;
	double a = p*q/(p+q);
	std::vector<double> P12 = P(exp1, exp2, xyz1, xyz2);
	std::vector<double> Q34 = P(exp3, exp4, xyz3, xyz4);
	std::vector<double> PQ({P12[0]-Q34[0],P12[1]-Q34[1],P12[2]-Q34[2]});
	double RPQ = sqrt(PQ[0]*PQ[0]+PQ[1]*PQ[1]+PQ[2]*PQ[2]);
	
	double sum = 0;
	for(int t = 0; t <= (L1[0]+L2[0]); t++){
		for(int u = 0; u <= (L1[1]+L2[1]); u++){
			for(int v = 0; v <= (L1[2]+L2[2]); v++){
				for(int tau = 0; tau <= (L3[0]+L4[0]); tau++){
					for(int nu = 0; nu <= (L3[1]+L4[1]); nu++){
						for(int phi = 0; phi <= (L3[2]+L4[2]); phi++){
							sum += 	E(L1[0],L2[0],t  ,exp1,exp2,xyz1[0]-xyz2[0])  *
					 			E(L1[1],L2[1],u  ,exp1,exp2,xyz1[1]-xyz2[1])  *
				       				E(L1[2],L2[2],v  ,exp1,exp2,xyz1[2]-xyz2[2])  *
				       				E(L3[0],L4[0],tau,exp3,exp4,xyz3[0]-xyz4[0])  *
				       				E(L3[1],L4[1],nu ,exp3,exp4,xyz3[1]-xyz4[1])  *
				       				E(L3[2],L4[2],phi,exp3,exp4,xyz3[2]-xyz4[2])  *
				       				R(0,t+tau,u+nu,v+phi,a,PQ[0],PQ[1],PQ[2],RPQ) *
								pow(-1, tau+nu+phi);
						}
					}
				} 
			}
		}
	}
	return 2*pow(M_PI, 2.5)/(p*q*sqrt(p+q)) * sum;	
}

double G(GF g1, GF g2, GF g3, GF g4){
        double sum = 0;
        for(int i = 0; i < g1.exps.size(); i++){
                for(int j = 0; j < g2.exps.size(); j++){
			for(int k = 0; k < g3.exps.size(); k++){
				for(int l = 0; l < g4.exps.size(); l++){
					sum +=  g1.N[i] * g2.N[j] * g3.N[k] * g4.N[l] * 
						g1.d[i] * g2.d[j] * g3.d[k] * g4.d[l] *
						Gp(g1.shell, g2.shell, g3.shell, g4.shell,
						   g1.exps[i], g2.exps[j], g3.exps[k], g4.exps[l],
						   g1.xyz, g2.xyz, g3.xyz, g4.xyz);
				}
			}
                }
        }
        return sum;
}
