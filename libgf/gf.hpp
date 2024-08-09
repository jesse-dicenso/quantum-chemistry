#ifndef GFHEADERDEF
#define GFHEADERDEF

#include "../libmath/miscmath.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class PGF{
	public:
		PGF(double exponent, std::vector<double> pos, std::vector<int> shl);
		PGF(double exponent, std::vector<double> pos, std::vector<int> shl, double Nrm);
		
		double exp;
		std::vector<double> xyz;
		std::vector<int> shell;
		double N;

		double getN();
		void setN(double nrm);

		friend bool operator== (const PGF &p1, const PGF &p2);

};

class CGF{
	public:
		CGF(std::vector<PGF> pgfC, std::vector<double> dC);
		
		int len;
		std::vector<PGF> pgf;
		std::vector<double> d;
		double N;
		
		friend bool operator== (const CGF &c1, const CGF &c2);
};

// Hermite Gaussian Expansion Coefficients
double E(int i, int j, int t, double a, double b, double QAB);

// components of coefficients
// for contractions of arbitrary angular momentum
// in Gaussian Product Theorem
double fk(int k, int y1, int y2, double PA, double PB);

// Product Constant
std::vector<double> K(PGF p1, PGF p2);

// New Center
std::vector<double> P(PGF p1, PGF p2);

#endif
