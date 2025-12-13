#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include "../libmol/mol.hpp"

#include <cmath>
#include <vector>

class grid{
	public:
		grid(const Molecule &mol, int num_radial = 100, int num_angular = 230, int becke_k = 3);

		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;
		std::vector<double> w;
		
		int num_gridpoints;
};

double becke_weight(const Molecule& mol, double x, double y, double z, int atom_me, int k = 3);
double becke_step(double mu, int k);

template <typename F, typename... Args>
double integrate_quad(const grid& g, F&& func, const Args&... args){
	double result = 0;
	for(int i = 0; i < g.num_gridpoints; i++){
		result += g.w[i] * func(g.x[i], g.y[i], g.z[i], args...);
	}
	return result;
}

// Only Lebedev grids of degree 230 available
extern const int 	lebedev_degree;
extern const double lebedev_x[230];
extern const double lebedev_y[230];
extern const double lebedev_z[230];
extern const double lebedev_w[230];

extern const double bragg_slater_radii[118];

#endif
