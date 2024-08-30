#ifndef ONEEHEADERDEF
#define ONEEHEADERDEF

#include "../libgf/gf.hpp"
#include "../libmath/linalg.hpp"

double Sp(std::vector<int> L1, std::vector<int> L2, double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2);
double S(GF g1, GF g2);

double Tp(std::vector<int> L1, std::vector<int> L2, double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2);
double T(GF g1, GF g2);

double R(int n, int t, int u, int v, double p, double XPC, double YPC, double ZPC, double RPC);

double Vp(std::vector<int> L1, std::vector<int> L2, double exp1, double exp2, std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyzN);
double V(GF g1, GF g2, std::vector<double> xyzN);

#endif
