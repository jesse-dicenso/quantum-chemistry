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

std::vector<double> density_gradient(double x, double y, double z, void* ctx){
	density_context* d_ctx = static_cast<density_context*>(ctx);
	const Molecule *mol = d_ctx->molecule;
	const Matrix *P = d_ctx->Pmatrix;

	std::vector<double> grad_rho = {0.0, 0.0, 0.0};
	std::vector<double> eval_grad_gf;
	double eval_gf;

	for(int i = 0; i < P->rows; i++){
			eval_grad_gf = mol->AOs[i].evaluate_gradient(x, y, z);
			for(int j = 0; j < P->cols; j++){
					eval_gf = P->matrix[i][j] * mol->AOs[j].evaluate(x, y, z);
					grad_rho[0] += eval_grad_gf[0] * eval_gf;
					grad_rho[1] += eval_grad_gf[1] * eval_gf;
					grad_rho[2] += eval_grad_gf[2] * eval_gf;
			}
	}
	grad_rho[0] *= 2;
	grad_rho[1] *= 2;
	grad_rho[2] *= 2;
	return grad_rho;
}

double integrate_density(const grid &g, const Molecule &mol, const Matrix &P){
	density_context ctx = {&mol, &P};
	return integrate_quad(g, density, &ctx);
}
