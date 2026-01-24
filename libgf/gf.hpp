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
		GF(const std::vector<double>& exponents, const std::vector<double>& coeffs, const std::vector<double>& pos, 
		   const std::vector<int>& shl, int atom_idx);
		
		std::vector<double> exps;
		std::vector<double> d;
		std::vector<double> xyz;
		std::vector<int> shell;
		
		int atom_index;
		
		std::vector<double> N;

		double evaluate(double x, double y, double z) const;
		std::vector<double> evaluate_gradient(double x, double y, double z) const;

		friend bool operator== (const GF &g1, const GF &g2);

};

// Hermite Gaussian Expansion Coefficients
double E(int i, int j, int t, double a, double b, double QAB);

// New Center
std::vector<double> P(double exp1, double exp2, const std::vector<double>& xyz1, const std::vector<double>& xyz2);

#endif
