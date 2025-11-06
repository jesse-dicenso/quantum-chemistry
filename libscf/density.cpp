#include "density.hpp"

struct R_context{
	const Molecule 	*molecule;
	const Matrix	*Pmatrix;
};

double R_density(double x, double y, double z, const Molecule &mol, const Matrix &P){
	assert(mol.NUPDOWN==0);
	double density = 0;
	double temp;

	// atomic
	if(mol.Natoms == 1){
		for(int i = 0; i < P.rows; i++){
			for(int j = i+1; j < P.cols; j++){
				density += 2 * P.matrix[i][j] * mol.AOs[i].evaluate(x, y, z) * mol.AOs[j].evaluate(x, y, z);
			}
			density += P.matrix[i][i] * pow(mol.AOs[i].evaluate(x, y, z), 2);
		}
		return density;
	}

	// molecular; use Becke partitioning
	else{
		return 0;
	}
}

double R_density_wrapper(double x, double y, double z, void* ctx){
	R_context* r_ctx = static_cast<R_context*>(ctx);
	return R_density(x, y, z, *r_ctx->molecule, *r_ctx->Pmatrix);
}

double integrate_R_density(const Molecule &mol, const Matrix &P, double b_s_r, int n){
	R_context ctx = {&mol, &P};
	double result = lebedev_gauss_chebyshev(R_density_wrapper, &ctx, mol.xyz[0][0], mol.xyz[0][1], mol.xyz[0][2], b_s_r, n);
	return result;
}
