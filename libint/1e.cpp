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

double R(int n, int t, int u, int v, double p, double XPC, double YPC, double ZPC, double RPC){
	assert((n>=0) && (t>=0) && (u>=0) && (v>=0));
	if((t+u+v)==0){
		return pow(-2*p, n) * boys(n, p*RPC*RPC);
	}
	else if(t!=0){
		if(t > 1){
			return (t-1)*R(n+1,t-2,u,v,p,XPC,YPC,ZPC,RPC) + XPC*R(n+1,t-1,u,v,p,XPC,YPC,ZPC,RPC);
		}
		else{
			return XPC*R(n+1,t-1,u,v,p,XPC,YPC,ZPC,RPC);
		}
	}
	else if(u!=0){
		if(u > 1){
			return (u-1)*R(n+1,t,u-2,v,p,XPC,YPC,ZPC,RPC) + YPC*R(n+1,t,u-1,v,p,XPC,YPC,ZPC,RPC);
		}
		else{
			return YPC*R(n+1,t,u-1,v,p,XPC,YPC,ZPC,RPC);
		}	
	}
	else{
		if(v > 1){
			return (v-1)*R(n+1,t,u,v-2,p,XPC,YPC,ZPC,RPC) + ZPC*R(n+1,t,u,v-1,p,XPC,YPC,ZPC,RPC);
		}
		else{
			return ZPC*R(n+1,t,u,v-1,p,XPC,YPC,ZPC,RPC);
		}
	}
}

double Vp(std::vector<int> L1, std::vector<int> L2, double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyzN){
	int i = L1[0];
	int j = L2[0];
	int k = L1[1];
	int l = L2[1];
	int m = L1[2];
	int n = L2[2];
	
	double p = exp1 + exp2;
	std::vector<double> P12 = P(exp1,exp2,xyz1,xyz2);
	std::vector<double> PC({P12[0]-xyzN[0],P12[1]-xyzN[1],P12[2]-xyzN[2]});
	double RPC = sqrt(PC[0]*PC[0]+PC[1]*PC[1]+PC[2]*PC[2]);	

	double sum = 0;
	for(int t = 0; t <= (i+j); t++){
		for(int u = 0; u <= (k+l); u++){
			for(int v = 0; v <= (m+n); v++){
				sum += E(i,j,t,exp1,exp2,xyz1[0]-xyz2[0]) *
				       E(k,l,u,exp1,exp2,xyz1[1]-xyz2[1]) *
				       E(m,n,v,exp1,exp2,xyz1[2]-xyz2[2]) *
				       R(0,t,u,v,p,PC[0],PC[1],PC[2],RPC);
			} 
		}
	}
	return (2*M_PI/p)*sum;
}

double V(GF g1, GF g2, std::vector<double> xyzN){
        double sum = 0;
        for(int i = 0; i < g1.exps.size(); i++){
                for(int j = 0; j < g2.exps.size(); j++){
			sum += g1.N[i] * g2.N[j] * g1.d[i] * g2.d[j] * Vp(g1.shell, g2.shell, g1.exps[i], g2.exps[j], g1.xyz, g2.xyz, xyzN);
                }
        }
        return sum;
}
