#include "density.hpp"

struct density_context{
	const Molecule 	*molecule;
	const Matrix	*Pmatrix;
};

double density(double x, double y, double z, void* ctx){
	density_context* d_ctx = static_cast<density_context*>(ctx);
	const Molecule *mol = d_ctx->molecule;
	const Matrix *P = d_ctx->Pmatrix;
	double rho = 0;
	double eval_gf;
	for(int i = 0; i < P->rows; i++){
			eval_gf = mol->AOs[i].evaluate(x, y, z);
			for(int j = i+1; j < P->cols; j++){
					rho += 2 * P->matrix[i][j] * eval_gf * mol->AOs[j].evaluate(x, y, z);
			}
			rho += P->matrix[i][i] * eval_gf * eval_gf;
	}
	return rho;
}

double integrate_density(const grid &g, const Molecule &mol, const Matrix &P){
	density_context ctx = {&mol, &P};
	return integrate_quad(g, density, &ctx);
}
