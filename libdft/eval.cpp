#include "eval.hpp"
#include "func.hpp"

// Evaluate E_XC, F_XC //////////////////////////////////////////

void LDA(XC* xc, XC_ret (*func)(XC*)){
	int num_gpts = xc->g->num_gridpoints;
	std::vector<double> phi_buf(xc->mol->AOs.size());
	zero_xc_data(xc);
	for(int gpt = 0; gpt < num_gpts; gpt++){
		eval_basis_funcs_per_gpt(xc, phi_buf, gpt);
		eval_density_per_gpt(xc, phi_buf);
		LDA_per_gpt(xc, phi_buf, func, gpt);
	}
}

// Helpers //////////////////////////////////////////////////////

void zero_xc_data(XC* xc){
	xc->E_XC = 0.0;
	if(xc->restricted){
		assert(xc->FXC!=nullptr);
		*(xc->FXC) = zero(xc->FXC->rows, xc->FXC->cols);
	}
	else{
		assert((xc->FXC_A!=nullptr) && (xc->FXC_B!=nullptr));
		*(xc->FXC_A) = zero(xc->FXC_A->rows, xc->FXC_A->cols);
		*(xc->FXC_B) = zero(xc->FXC_B->rows, xc->FXC_B->cols);
	}
}

void eval_basis_funcs_per_gpt(XC* xc, std::vector<double>& phi_buf, int gpix){
	for(int j = 0; j < xc->mol->AOs.size(); j++){
		phi_buf[j] = xc->mol->AOs[j].evaluate(xc->g->x[gpix], xc->g->y[gpix], xc->g->z[gpix]);
	}	
}

void eval_density_per_gpt(XC* xc, const std::vector<double>& phi_buf){
	if(xc->restricted){xc->rho = density(phi_buf, *xc->P);}
	else{
		xc->rho_a = density(phi_buf, *xc->P_A);
		xc->rho_b = density(phi_buf, *xc->P_B);
		xc->rho = xc->rho_a + xc->rho_b;
	}
}

void LDA_per_gpt(XC* xc, const std::vector<double>& phi_buf, XC_ret (*func)(XC*), int gpix){
	Matrix* fxc[2];
	fxc[0] = (xc->restricted ? xc->FXC : xc->FXC_A);
	fxc[1] = (xc->restricted ? nullptr : xc->FXC_B);

	int dim   = (xc->restricted ? xc->FXC->rows : xc->FXC_A->rows);
	int spins = (xc->restricted ? 1 : 2);

	if(xc->rho < 1e-20){return;}
	XC_ret ret = func(xc);
	xc->E_XC += xc->g->w[gpix] * ret.e_XC;
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			fxc[s]->matrix[mu][mu] += xc->g->w[gpix] * phi_buf[mu] * ret.v_XC[s] * phi_buf[mu];
			for(int nu = 0; nu < mu; nu++){
				fxc[s]->matrix[mu][nu] += xc->g->w[gpix] * phi_buf[mu] * ret.v_XC[s] * phi_buf[nu];
				fxc[s]->matrix[nu][mu] = fxc[s]->matrix[mu][nu];
			}
		}
	}
}
