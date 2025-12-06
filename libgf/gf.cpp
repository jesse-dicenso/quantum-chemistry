#include "gf.hpp"

// class GF

GF::GF(std::vector<double> exponents, std::vector<double> coeffs, std::vector<double> pos, std::vector<int> shl, int atom_idx){
	assert(exponents.size()==coeffs.size());
	assert(pos.size()==3);
	assert(shl.size()==3);
	assert((shl[0] >= 0) && (shl[1] >= 0) && (shl[2] >=0));
	
	exps = exponents;
	d = coeffs;
	xyz = pos;
	shell = shl;
	atom_index = atom_idx;

	int l = shell[0];
	int m = shell[1];
	int n = shell[2];
	int L = l+m+n;
	
	// Normalize
	std::vector<double> nrms(exps.size());
	for(int i = 0; i < exps.size(); i++){
		nrms[i] = pow((2/M_PI),0.75) * intpow(2,L) * pow(exps[i],(2*L+3)/4.0) / sqrt(dfact(2*l-1)*dfact(2*m-1)*dfact(2*n-1));
	}
	
	N = nrms;
	
	double cN = 0;
	for(int i = 0; i < exps.size(); i++){
		for(int j = 0; j < exps.size(); j++){
			cN += N[i]*N[j]*d[i]*d[j] / pow((exps[i] + exps[j]), (1.5 + L));
		}
	}
	
	cN = 1 / sqrt(pow(M_PI, 1.5) * dfact(2*l-1) * dfact(2*m-1) * dfact(2*n-1) * cN / intpow(2, L));	
	
	for(int i = 0; i < exps.size(); i++){
		d[i] *= cN;
	}
}

double GF::evaluate(double x, double y, double z) const{
	double sum = 0;
	double r2 = (x-xyz[0])*(x-xyz[0]) + (y-xyz[1])*(y-xyz[1]) + (z-xyz[2])*(z-xyz[2]);
	for(int i = 0; i < exps.size(); i++){
		sum += N[i] * d[i] * exp(-exps[i]*r2);
	}
	return sum * intpow(x-xyz[0], shell[0]) * intpow(y-xyz[1], shell[1]) * intpow(z-xyz[2], shell[2]);
}
		
std::vector<double> GF::evaluate_gradient(double x, double y, double z) const{
	std::vector<double> gradient = {0.0, 0.0, 0.0};
	double sum = 0;
	double r2 = (x-xyz[0])*(x-xyz[0]) + (y-xyz[1])*(y-xyz[1]) + (z-xyz[2])*(z-xyz[2]);
	double factor_x = intpow(x-xyz[0], shell[0]);
	double factor_y = intpow(y-xyz[1], shell[1]);
	double factor_z = intpow(z-xyz[2], shell[2]);
	for(int i = 0; i < exps.size(); i++){
		sum -= N[i] * d[i] * exps[i] * exp(-exps[i]*r2);
	}
	sum *= 2 * factor_x * factor_y * factor_z;
	gradient[0] = sum * (x-xyz[0]);
	gradient[1] = sum * (y-xyz[1]);
	gradient[2] = sum * (z-xyz[2]);

	bool isl = (shell[0]!=0);
	bool ism = (shell[1]!=0);
	bool isn = (shell[2]!=0);
	if(isl || ism || isn){
		sum = 0;	
		for(int i = 0; i < exps.size(); i++){
			sum += N[i] * d[i] * exp(-exps[i]*r2);
		}
		if(isl) {gradient[0] += shell[0] * intpow(x-xyz[0], shell[0]-1) * factor_y * factor_z * sum;}
		if(ism) {gradient[1] += shell[1] * intpow(y-xyz[1], shell[1]-1) * factor_z * factor_x * sum;}
		if(isn) {gradient[2] += shell[2] * intpow(z-xyz[2], shell[2]-1) * factor_x * factor_y * sum;}
	}
	return gradient;
}

bool operator== (const GF &g1, const GF &g2){
	if((g1.exps==g2.exps) && (g1.d==g2.d) && (g1.shell==g2.shell) && (g1.atom_index==g2.atom_index)){
		return true;
	}
	else{
		return false;
	}
}

// misc functions related to GF

double E(int i, int j, int t, double a, double b, double QAB){
	assert((i>=0) && (j>=0));
	double p = a + b;
	if((t < 0) || (t > (i+j))){
		return 0;
	}
	else if((t==0) && ((i+j)==0)){
		return exp(-a*b*QAB*QAB/p);
	}
	else if(j==0){
		return E(i-1, j, t-1, a, b, QAB)/(2*p) - (b*QAB/p)*E(i-1, j, t, a, b, QAB) + (t+1)*E(i-1, j, t+1, a, b, QAB);
	}
	else{
		return E(i, j-1, t-1, a, b, QAB)/(2*p) + (a*QAB/p)*E(i, j-1, t, a, b, QAB) + (t+1)*E(i, j-1, t+1, a, b, QAB);
	}
}

std::vector<double> P(double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2){
	std::vector<double> vec(3);
	for(int i = 0; i < 3; i++){
		vec[i] = (exp1*xyz1[i]+exp2*xyz2[i])/(exp1+exp2);
	}
	return vec;
}
