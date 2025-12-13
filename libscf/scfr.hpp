#ifndef SCFRHEADERDEF
#define SCFRHEADERDEF

#include "scfgen.hpp"
#include "../libdft/xc.hpp"

Matrix R_density_matrix(const Matrix& C, int N);
double R_E0(XC_inp* xc_inp, const Matrix& Hcore, const Matrix& F, const Matrix& J);

void R_FPI (const Matrix& s, const Matrix& hcore, const Matrix& x, XC_inp* xc_inp, Matrix* f, Matrix* fo, 
			Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i);
void R_DIIS(const Matrix& s, const Matrix& hcore, const Matrix& x, XC_inp* xc_inp, Matrix* f, Matrix* fo, 
			Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, std::vector<Matrix>& SPf, 
			std::vector<Matrix>& SPe, int sps, int* icd);

#endif
