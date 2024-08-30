#ifndef GFHEADERDEF
#define GFHEADERDEF

#include "../libmath/miscmath.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class GF{
	public:
		GF(std::vector<double> exponents, std::vector<double> coeffs, std::vector<double> pos, std::vector<int> shl);
		
		std::vector<double> exps;
		std::vector<double> d;
		std::vector<double> xyz;
		std::vector<int> shell;
		
		std::vector<double> N;

		void setN(double nrm);

		friend bool operator== (const GF &g1, const GF &g2);

};

// Hermite Gaussian Expansion Coefficients
double E(int i, int j, int t, double a, double b, double QAB);

// New Center
std::vector<double> P(double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2);

#endif
