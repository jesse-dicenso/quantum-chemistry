#ifndef TWOEHEADERDEF
#define TWOEHEADERDEF

#include "1e.hpp"

class ERI{
	public:
		ERI(int p, int q, int r, int s);
		
		int a;
		int b;
		int c;
		int d;
		double eri;

		bool isEqual(ERI e2);
};

double Gp(std::vector<int> L1, std::vector<int> L2, std::vector<int> L3, std::vector<int> L4, 
	  double exp1, double exp2, double exp3, double exp4,
	  std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyz3, std::vector<double> xyz4);
double G(GF g1, GF g2, GF g3, GF g4);

std::vector<std::vector<std::vector<std::vector<double>>>> ERIs(std::vector<GF> phis);

#endif
