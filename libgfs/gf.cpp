#include "gf.hpp"

// class PGF

PGF::PGF(double exponent, std::vector<double> pos, int Ll, int Lm, int Ln){
	assert((Ll >= 0) && (Lm >= 0) && (Ln >=0));
	assert(pos.size()==3);
	exp = exponent;
	xyz = pos;
	l = Ll;
	m = Lm;
	n = Ln;

	if((l==0) && (m==0) && (n==0)){
		N = pow(2.0 * exp / M_PI, 0.75);
	}
	else{
		double dfl = 1.0;
		double dfm = 1.0;
		double dfn = 1.0;
		if(l > 0){
			dfl = dfact(2*l-1);
		}
		if(m > 0){
			dfm = dfact(2*m-1);
		}
		if(n > 0){
			dfn = dfact(2*n-1);
		}
		N = pow(2.0 * exp / M_PI, 0.75) * pow(pow(4*exp, l+m+n)/(dfl*dfm*dfn), 0.5);
	}
}

int PGF::getl(){
	return l;
}

int PGF::getm(){
	return m;
}

int PGF::getn(){
	return n;
}

double PGF::getN(){
	return N;
}

double PGF::getexp(){
	return exp;
}

void PGF::printpgf(){
	std::cout << "Start PGF\n";
	std::cout << "   Exponent: " << exp << '\n';
	std::cout << "   Center (x, y, z): " << std::setprecision(15) << xyz[0] << ", " << std::setprecision(15) << xyz[1] << ", " << std::setprecision(15) << xyz[2] << '\n';
	std::cout << "   Angular momenta (l, m, n): " << l << ", " << m << ", " << n << '\n';
	std::cout << "   Normalization constant: " << N << '\n';
	std::cout << "End PGF\n";
}

bool operator== (const PGF &p1, const PGF &p2){
	if((p1.exp==p2.exp) && (p1.xyz==p2.xyz) && (p1.l==p2.l) && (p1.m==p2.m) && (p1.n==p2.n)){
		return true;
	}
	else{
		return false;
	}
}

// class CGF

CGF::CGF(int lenC, std::vector<PGF> pgfC, std::vector<double> dC){
	assert((lenC == pgfC.size()) && (lenC == dC.size()));
	// Center and angular momenta of each PGF must be the same!
	for(int h = 1; h < lenC; h++){
		assert((pgfC[0].xyz[0]==pgfC[h].xyz[0]) && (pgfC[0].xyz[1]==pgfC[h].xyz[1]) && (pgfC[0].xyz[2]==pgfC[h].xyz[2]));
		assert((pgfC[0].getl()==pgfC[h].getl()) && (pgfC[0].getm()==pgfC[h].getm()) && (pgfC[0].getn()==pgfC[h].getn()));
	}
		
	len = lenC;
	pgf = pgfC;
	d = dC;
	if((pgf[0].getl()==0) && (pgf[0].getm()==0) && (pgf[0].getn()==0)){
		double sum = 0.0;
		for(int i = 0; i < len; i++){
			for(int j = 0; j < len; j++){
				sum += d[i]*d[j] / pow((pgf[i].getexp() + pgf[j].getexp()), 1.5);
			}
		}
		N = pow(M_PI, -0.75) / pow(sum, 0.5);
	}
	else{
		double sum = 0.0;
		double temp;
		double lt = pgf[0].getl();
		double mt = pgf[0].getm();
		double nt = pgf[0].getn();
		double dfl = 1.0;
		double dfm = 1.0;
		double dfn = 1.0;
		if(lt > 0){
			dfl = dfact(2*lt-1);
		}
		if(mt > 0){
			dfm = dfact(2*mt-1);
		}
		if(nt > 0){
			dfn = dfact(2*nt-1);
		}
		for(int i = 0; i < len; i++){
			for(int j = 0; j < len; j++){
				sum += d[i]*d[j] / pow((pgf[i].getexp() + pgf[j].getexp()), (1.5 + lt + mt + nt));
			}
		}
		temp = pow(pow(2, (lt + mt + nt)) / (pow(M_PI, 1.5) * dfl * dfm * dfn), 0.5);
		N = temp / pow(sum, 0.5);
	}
}

int CGF::getlen(){
	return len;
}

std::vector<PGF> CGF::getpgf(){
	return pgf;
}

std::vector<double> CGF::getd(){
	return d;
}

double CGF::getN(){
	return N;
}

void CGF::printcgf(){
	std::cout << "Start CGF\n";
	std::cout << "Contraction Length: " << len << '\n';
	std::cout << "Normalization constant: " << std::setprecision(15) << N << '\n';
	std::cout << "Center (x, y, z): " << std::setprecision(15) << pgf[0].xyz[0] << ", " << std::setprecision(15) << pgf[0].xyz[1] << ", "  << std::setprecision(15) << pgf[0].xyz[2] << '\n';
	std::cout << "Angular momenta (l, m, n): " << pgf[0].getl() << ", " << pgf[0].getm() << ", " << pgf[0].getn() << '\n';
	for(int i = 0; i < len; i++){
		std::cout << "   PGF " << i << '\n';
		std::cout << "      d = " << d[i] << '\n';
		std::cout << "      Exponent: " << pgf[i].getexp() << '\n';
	}
	std::cout << "End CGF\n";
}

bool operator== (const CGF &c1, const CGF &c2){
	if((c1.len == c2.len) && (c1.N == c2.N) && (c1.d == c2.d) && (c1.pgf == c2.pgf)){
		return true;
	}
	else{
		return false;
	}
}

// misc functions related to PGF/CGF

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

double K(PGF p1, PGF p2){
	std::vector<double> AB(3);
	for(int i = 0; i < 3; i++){
		AB[i] = p1.xyz[i] - p2.xyz[i];
	}
	return exp(-(p1.getexp()*p2.getexp()/(p1.getexp()+p2.getexp()))*dot(AB,AB));
}

std::vector<double> P(PGF p1, PGF p2){
	std::vector<double> vec(3);
	double a = p1.getexp();
	double b = p2.getexp();
	for(int i = 0; i < 3; i++){
		vec[i] = (a*p1.xyz[i]+b*p2.xyz[i])/(a+b);
	}
	return vec;
}
