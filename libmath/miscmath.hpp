#ifndef MISCMATHHEADERDEF
#define MISCMATHHEADERDEF

#include <cassert>
#include <cmath>
#include <vector>

double intpow(double x, int n);
double long fact(double long n);
double long dfact(double long n);
double dot(const std::vector<double>& a, const std::vector<double>& b);
double boys(int n, double x);

#endif
