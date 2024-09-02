#ifndef SCFRHFHEADERDEF
#define SCFRHFHEADERDEF

#include "../libint/2e.hpp"
#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

Matrix R_density_matrix(Matrix C, int N);
Matrix R_F(Matrix Hcore, Matrix P, std::vector<std::vector<std::vector<std::vector<double>>>> eris);
double R_E0(Matrix P, Matrix Hcore, Matrix F);

void R_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N);
void R_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps, int* icd);

#endif
