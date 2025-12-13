#ifndef SCFUHEADERDEF
#define SCFUHEADERDEF

#include "scfgen.hpp"
#include "../libdft/xc.hpp"

Matrix UR_density_matrix(const Matrix& C, int N);
double UR_E0(XC_inp* xc_inp, const Matrix& Hcore, const Matrix& Fa, const Matrix& Fb, const Matrix& J);

void UR_FPI (const Matrix& s, const Matrix& hcore, const Matrix& x, XC_inp* xc_inp, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, 
			 Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i);
void UR_DIIS(const Matrix& s, const Matrix& hcore, const Matrix& x, XC_inp* xc_inp, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, 
			 Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, 
			 std::vector<Matrix>& SPfa, std::vector<Matrix>& SPfb, std::vector<Matrix>& SPea, std::vector<Matrix>& SPeb, int sps, int* icd);

#endif
