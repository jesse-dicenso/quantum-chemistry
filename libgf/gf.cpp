#include "gf.hpp"

// class PGF

PGF::PGF(double exponent, std::vector<double> pos, std::vector<int> shl){
	assert(pos.size()==3);
	assert(shl.size()==3);
	assert((shl[0] >= 0) && (shl[1] >= 0) && (shl[2] >=0));
	exp = exponent;
	xyz = pos;
	shell = shl;
		
	N = pow((2/M_PI),0.75) * pow(2,(shell[0]+shell[1]+shell[2])) * pow(exp,(2*(shell[0]+shell[1]+shell[2])+3)/4.0) / sqrt(dfact(2*shell[0]-1)*dfact(2*shell[1]-1)*dfact(2*shell[2]-1));
}

PGF::PGF(double exponent, std::vector<double> pos, std::vector<int> shl, double Nrm){
	assert(pos.size()==3);
	assert(shl.size()==3);
	assert((shl[0] >= 0) && (shl[1] >= 0) && (shl[2] >=0));
	exp = exponent;
	xyz = pos;
	shell = shl;
	N = Nrm;
}

bool operator== (const PGF &p1, const PGF &p2){
	if((p1.exp==p2.exp) && (p1.xyz==p2.xyz) && (p1.shell==p2.shell)){
		return true;
	}
	else{
		return false;
	}
}

// class CGF

CGF::CGF(std::vector<PGF> pgfC, std::vector<double> dC){
	assert(pgfC.size() == dC.size());
	// Center and angular momenta of each PGF must be the same!
	for(int h = 1; h < len; h++){
		for(int f = 0; f < 3; f++){	
			assert(pgfC[0].xyz[f]==pgfC[h].xyz[f]);
			assert(pgfC[0].shell[f]==pgfC[h].shell[f]);
		}
	}

	len = pgfC.size();
	pgf = pgfC;
	d = dC;

	double sum = 0;
	for(int i = 0; i < len; i++){
		for(int j = 0; j < len; j++){
			sum += d[i]*d[j] / pow((pgf[i].exp + pgf[j].exp), (1.5 + pgf[0].shell[0]+pgf[0].shell[1]+pgf[0].shell[2]));
		}
	}
	N = 1 / sqrt(pow(M_PI, 1.5) * dfact(2*pgf[0].shell[0]-1) 
				    * dfact(2*pgf[0].shell[1]-1)
				    * dfact(2*pgf[0].shell[2]-1) 
				    * sum / pow(2, (pgf[0].shell[0]+pgf[0].shell[1]+pgf[0].shell[2])));
	for(int j = 0; j < len; j++){
		d[j] *= N;
	}
}

bool operator== (const CGF &c1, const CGF &c2){
	if((c1.len == c2.len) && (c1.d == c2.d) && (c1.pgf == c2.pgf)){
		return true;
	}
	else{
		return false;
	}
}

// misc functions related to PGF/CGF

double E(int i, int j, int t, double a, double b, double QAB){
	assert((i>=0) && (j>=0));
	double p = a + b;
	if((t < 0) || (t > (i+j))){
		return 0;
	}
	else if((t==0) || ((i+j)==0)){
		return exp(-a*b*QAB*QAB/p);
	}
	else if(j==0){
		return E(i-1, j, t-1, a, b, QAB)/(2*p) - (b*QAB/p)*E(i-1, j, t, a, b, QAB) + (t+1)*E(i-1, j, t+1, a, b, QAB);
	}
	else{
		return E(i, j-1, t-1, a, b, QAB)/(2*p) + (a*QAB/p)*E(i, j-1, t, a, b, QAB) + (t+1)*E(i, j-1, t+1, a, b, QAB);
	}
}

double fk(int k, int y1, int y2, double PA, double PB){
        double sum = 0;
        for(int i = 0; i <= y1; i++){
                for(int j = 0; j <= y2; j++){
                        if((i+j)==k){
                                sum += pow(PA,(y1-i))*pow(PB,(y2-j))*binomial(y1,i)*binomial(y2,j);
                        }
                }
        }
        return sum;
}

std::vector<double> K(PGF p1, PGF p2){
	std::vector<double> vec(3);
	double mu = p1.exp*p2.exp/(p1.exp+p2.exp);
	for(int i = 0; i < 3; i++){
		vec[i] = exp(-mu*(p1.xyz[i]-p2.xyz[i])*(p1.xyz[i]-p2.xyz[i]));
	}
	return vec;
}

std::vector<double> P(PGF p1, PGF p2){
	std::vector<double> vec(3);
	for(int i = 0; i < 3; i++){
		vec[i] = (p1.exp*p1.xyz[i]+p2.exp*p2.xyz[i])/(p1.exp+p2.exp);
	}
	return vec;
}
