#ifndef SCFUHFHEADERDEF
#define SCFUHFHEADERDEF

#include "../libint/2e.hpp"
#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

Matrix UR_density_matrix(Matrix C, int N);
Matrix UR_F(Matrix Hcore, Matrix PT, Matrix Ps, std::vector<std::vector<std::vector<std::vector<double>>>> g);
double UR_E0(Matrix PT, Matrix Pa, Matrix Pb, Matrix Hcore, Matrix Fa, Matrix Fb);

void UR_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb);
void UR_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, Matrix* SPfa, Matrix* SPfb, Matrix* SPea, Matrix* SPeb, int sps, int* icd);

#endif
