#ifndef SCFUHEADERDEF
#define SCFUHEADERDEF

#include "scfgen.hpp"

Matrix UR_density_matrix(const Matrix& C, int N);
Matrix UR_F(const Matrix& Hcore, const Matrix& PT, const Matrix& Ps, const std::vector<std::vector<std::vector<std::vector<double>>>>& g);
double UR_E0(const Matrix& PT, const Matrix& Pa, const Matrix& Pb, const Matrix& Hcore, const Matrix& Fa, const Matrix& Fb);

void UR_FPI (const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, 
			 const Matrix& x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, 
			 Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i);
void UR_DIIS(const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, 
			 const Matrix& x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, 
			 Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, 
			 int i, std::vector<Matrix>& SPfa, std::vector<Matrix>& SPfb, std::vector<Matrix>& SPea, std::vector<Matrix>& SPeb, 
			 int sps, int* icd);

#endif
