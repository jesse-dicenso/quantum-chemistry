#ifndef BUILDMATRIXHEADERDEF
#define BUILDMATRIXHEADERDEF

#include "2e.hpp"
#include "../libmath/linalg.hpp"

Matrix overlap(std::vector<GF> phis);
Matrix kinetic(std::vector<GF> phis);
Matrix nuclear(std::vector<GF> phis, std::vector<int> Zvals, std::vector<std::vector<double>> xyzN);

#endif
