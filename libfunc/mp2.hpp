#ifndef MP2HEADERDEF
#define MP2HEADERDEF

#include "../libmath/linalg.hpp"

#include <vector>

double MP2_ENERGY(const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, const Matrix& C, const Matrix& E);

#endif
