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
		PGF(double exponent, std::vector<double> pos, int Ll, int Lm, int Ln);
		
		std::vector<double> xyz;

		double getexp();
		int getl();
		int getm();
		int getn();
		double getN();

		void printpgf();
		friend bool operator== (const PGF &p1, const PGF &p2);

	private:
		double exp;
		int l;
		int m;
		int n;
		double N;
};

class CGF{
	public:
		CGF(int lenC, std::vector<PGF> pgfC, std::vector<double> dC);
		
		int getlen();
		std::vector<PGF> getpgf();
		std::vector<double> getd();
		double getN();
		void printcgf();
		friend bool operator== (const CGF &c1, const CGF &c2);

	private:
		int len;
		std::vector<PGF> pgf;
		std::vector<double> d;
		double N;
};

// Functions for Gaussian Product Theorem

// fk: returns components of coefficients
// for contractions of arbitrary angular momentum
// in Gaussian Product Theorem
double fk(int k, int y1, int y2, double PA, double PB);

// Product Constant
double K(PGF p1, PGF p2);

// New Center
std::vector<double> P(PGF p1, PGF p2);

/*
double overlap(PGF A, PGF B);
double overlap(CGF A, CGF B);
*/

#endif
