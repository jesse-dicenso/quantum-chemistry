#ifndef SNXHELPERHEADERDEF
#define SNXHELPERHEADERDEF

#include "../libgf/gf.hpp"
#include "../libmath/linalg.hpp"
#include "../libmath/miscmath.hpp"

#include <cmath>
#include <vector>

class XC;

// Build 3 index tensor A_{\nu \lambda g} for batch of gridpoints
void SNX_A(const XC& xc, Tensor3& A, int g_start, int g_end);

// Contract A and F to give G^T, assumes g is slow idx of A
Matrix contract_A_F(const Tensor3& A, const Matrix& F);

// Misc functions for integral screening
double V_screen();
double E_rstu(int r, int s, int t, int u, double a, double b, double QAB);
double Theta(int k, int l, double a, double b, double QAB);

#endif
