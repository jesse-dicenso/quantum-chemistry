#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "../libint/2e.hpp"
#include "../libmath/linalg.hpp"
#include "../libmol/mol.hpp"

Matrix overlap(std::vector<GF> phis);
Matrix kinetic(std::vector<GF> phis);
Matrix nuclear(std::vector<GF> phis, std::vector<int> Zvals, std::vector<std::vector<double>> xyzN);
double nucrepl(std::vector<int> Z, std::vector<std::vector<double>> xyzN);

Matrix R_density_matrix(Matrix C, int N);
Matrix R_F(Matrix Hcore, Matrix P, std::vector<std::vector<std::vector<std::vector<double>>>> eris);
double R_E0(Matrix P, Matrix Hcore, Matrix F);

void R_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N);
void R_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps);

Matrix UR_density_matrix(Matrix C, int N);
Matrix UR_F(Matrix Hcore, Matrix PT, Matrix Ps, std::vector<std::vector<std::vector<std::vector<double>>>> g);
double UR_E0(Matrix PT, Matrix Pa, Matrix Pb, Matrix Hcore, Matrix Fa, Matrix Fb);

void UR_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int N);
void UR_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* foa, Matrix* fob, Matrix* ea, Matrix* eb, Matrix* coa, Matrix* cob, Matrix* ca, Matrix* cb, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps);

#endif
