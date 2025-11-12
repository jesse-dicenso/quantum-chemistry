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

extern const double PI_2p5;
double Gp(const std::vector<int>& L1, const std::vector<int>& L2, const std::vector<int>& L3, const std::vector<int>& L4, double exp1, double exp2, double exp3, double exp4, const std::vector<double>& xyz1, const std::vector<double>& xyz2, const std::vector<double>& xyz3, const std::vector<double>& xyz4);
double G(const GF& g1, const GF& g2, const GF& g3, const GF& g4);

std::vector<std::vector<std::vector<std::vector<double>>>> ERIs(const std::vector<GF>& phis);

#endif
