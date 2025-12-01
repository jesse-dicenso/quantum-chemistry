#include "density.hpp"

struct R_density_context{
	const Molecule 	*molecule;
	const Matrix	*Pmatrix;
};

double R_density(double x, double y, double z, void* ctx){
	R_density_context* r_ctx = static_cast<R_density_context*>(ctx);
	const Molecule *mol = r_ctx->molecule;
	assert(mol->NUPDOWN==0);
	const Matrix *P = r_ctx->Pmatrix;
	double density = 0;
	double eval_gf;
	for(int i = 0; i < P->rows; i++){
			eval_gf = mol->AOs[i].evaluate(x, y, z);
			for(int j = i+1; j < P->cols; j++){
					density += 2 * P->matrix[i][j] * eval_gf * mol->AOs[j].evaluate(x, y, z);
			}
			density += P->matrix[i][i] * eval_gf * eval_gf;
	}
	return density;
}

double integrate_R_density(const grid &g, const Molecule &mol, const Matrix &P){
	R_density_context ctx = {&mol, &P};
	return integrate_quad(g, R_density, &ctx);
}
