#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "../libint/2e.hpp"
#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

Matrix density_matrix(Matrix C, int N);
double E0(Matrix P, Matrix Hcore, Matrix F);

void FPI(Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* g, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N);
void DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* g, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps);

#endif
