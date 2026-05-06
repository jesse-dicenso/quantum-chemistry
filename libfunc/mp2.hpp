#ifndef MP2HEADERDEF
#define MP2HEADERDEF

#include "../libmath/linalg.hpp"

#include <vector>

double MP2_ENERGY(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, // AO ERIs
        const Matrix& C, // MO coefficient matrix, K * K
        const Matrix& e, // Orbital energies (HF eigenvalues), K * K
        size_t N,        // # electrons
        size_t K         // # basis functions
    );

#endif
