#ifndef SCFRHEADERDEF
#define SCFRHEADERDEF

#include "scfgen.hpp"

Matrix R_density_matrix(const Matrix& C, int N);
double R_E0(const Matrix& P, const Matrix& Hcore, const Matrix& F);

void R_FPI (const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, 
			const Matrix& x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, 
			int N, int i);
void R_DIIS(const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, 
			const Matrix& x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, 
			int N, int i, std::vector<Matrix>& SPf, std::vector<Matrix>& SPe, int sps, int* icd);

#endif
