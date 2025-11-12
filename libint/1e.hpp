#ifndef ONEEHEADERDEF
#define ONEEHEADERDEF

#include "../libgf/gf.hpp"
#include "../libmath/linalg.hpp"

double Sp(const std::vector<int>& L1, const std::vector<int>& L2, double exp1, double exp2, const std::vector<double>& xyz1, const std::vector<double>& xyz2);
double S(GF g1, GF g2);

double Tp(const std::vector<int>& L1, const std::vector<int>& L2, double exp1, double exp2, const std::vector<double>& xyz1, const std::vector<double>& xyz2);
double T(GF g1, GF g2);

double R(int n, int t, int u, int v, double p, double XPC, double YPC, double ZPC, double RPC);

double Vp(const std::vector<int>& L1, const std::vector<int>& L2, double exp1, double exp2, const std::vector<double>& xyz1, const std::vector<double>& xyz2, const std::vector<double>& xyzN);
double V(GF g1, GF g2, const std::vector<double>& xyzN);

#endif
