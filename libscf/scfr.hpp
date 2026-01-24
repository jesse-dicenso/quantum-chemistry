#ifndef SCFRHEADERDEF
#define SCFRHEADERDEF

#include "scfgen.hpp"
#include "../libdft/func.hpp"

Matrix R_density_matrix(const Matrix& C, int N);

void R_FPI (const Matrix& s, const Matrix& hcore, const Matrix& x, XC* xc, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, 
			double* Eo, double* err, int N, int i);
void R_DIIS(const Matrix& s, const Matrix& hcore, const Matrix& x, XC* xc, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, 
			double* Eo, double* err, int N, int i, std::vector<Matrix>& SPf, std::vector<Matrix>& SPe, int sps, int* icd);

#endif
