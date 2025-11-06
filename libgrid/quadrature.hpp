#ifndef QUADRATUREHEADERDEF
#define QUADRATUREHEADERDEF

#include<cmath>

double lebedev_gauss_chebyshev(double (*func)(double, double, double, void*), void* ctx, double c_x, double c_y, double c_z, double bragg_slater_r, int gauss_chebyshev_n);

// Lebedev points and weights for n = 25 (degree 230)
// Weights already normalized by 4PI
extern const int    lebedev_degree;
extern const double lebedev_x[230];
extern const double lebedev_y[230];
extern const double lebedev_z[230];
extern const double lebedev_w[230];

#endif
