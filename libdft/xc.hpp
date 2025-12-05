#ifndef XCHEADERDEF
#define XCHEADERDEF

#include "density.hpp"

// Exchange only
Matrix R_HF_X  (const Matrix& P , const std::vector<std::vector<std::vector<std::vector<double>>>>& g);
Matrix U_HF_X_s(const Matrix& Ps, const std::vector<std::vector<std::vector<std::vector<double>>>>& g);

// Correlation only
//Matrix R_Slater_X();
//Matrix U_Slater_X_s();

// Exchange and correlation

#endif
