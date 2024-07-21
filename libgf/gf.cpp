#include "gf.hpp"
#include "../libmath/miscmath.hpp" // Need for dfact(n)
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

PGF::PGF(double exponent, double Rx, double Ry, double Rz, int Ll, int Lm, int Ln){
	assert((Ll >= 0) && (Lm >= 0) && (Ln >=0));
	exp = exponent;
	x = Rx;
	y = Ry;
	z = Rz;
	l = Ll;
	m = Lm;
	n = Ln;

	if((l==0) && (m==0) && (n==0)){
		N = pow(2.0 * exp / M_PI, 0.75);
	}
	else{
		N = pow(2.0 * exp / M_PI, 0.75) * pow(pow(4*exp, l+m+n)/(dfact(2*l-1)*dfact(2*m-1)*dfact(2*n-1)), 0.5);
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

void PGF::printpgf(){
	cout << "*********************************\n";
	cout << "*  Primitive Gaussian Function  *\n";
	cout << "*********************************\n";
	cout << "Exponent:\n";
	cout << exp << '\n';
	cout << "Center (x, y, z):\n";
	cout << setprecision(15) << x << ", " << setprecision(15) << y << ", " << setprecision(15) << z << '\n';
	cout << "Angular momenta (l, m, n):\n";
	cout << l << ", " << m << ", " << n << '\n';
	cout << "Normalization constant:";
	cout << N << '\n';
	cout << "end PGF\n";
}

CGF::CGF(int lenC, vector<PGF> pgfC, vector<double> dC]){
	assert((lenC == pgfC.size()) && (lenC == dC.size()));
	len = lenC;
	pgf = pgfC;
	d = dC;

	// Angular momenta for each PGF must be the same!
	if((pgf[0].getl()==0) && (pgf[0].getm()==0) && (pgf[0].getn()==0)){
		double sum = 0.0;
		for(int i = 0; i < len; i++){
			for(int j = 0; j < len; j++){
				sum += d[i]*d[j] / pow((pgf[i].getexp() + pgf[j].getexp()), 1.5);
			}
		}
		N = pow(M_PI, 0.75) / pow(sum, 0.5);
	}
	else{
		double sum = 0.0;
		double temp;
		double lt = pgf[0].getl();
		double mt = pgf[0].getm();
		double nt = pgf[0].getn();
		for(int i = 0; i < len; i++){
			for(int j = 0; j < len; j++){
				sum += d[i]*d[j] / pow((pgf[i].getexp() + pgf[j].getexp()), (1.5 + lt + mt + nt));
			}
		}
		temp = pow(pow(2, (lt + mt + nt)) / (pow(M_PI, 1.5) * dfact(2*lt-1) * dfact(2*mt-1) * dfact(2*nt-1)), 0.5);
		N = temp / pow(sum, 0.5);
	}
}
