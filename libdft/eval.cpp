#include "eval.hpp"
#include "func.hpp"

// Evaluate E_XC, F_XC //////////////////////////////////////////

void LDA(XC* xc, LDA_ret (*func)(XC*)){
	std::vector<double> phi_buf(xc->mol->AOs.size());
	zero_xc_data(xc);
	for(int gpt = 0; gpt < xc->g->num_gridpoints; gpt++){
		eval_bfs_per_gpt(xc, phi_buf, gpt);
		eval_density_per_gpt(xc, phi_buf);
		LDA_per_gpt(xc, func, phi_buf, gpt);
	}
}

void GGA(XC* xc, GGA_ret (*func)(XC*)){
	std::vector<double> phi_buf(xc->mol->AOs.size());
	std::vector<double> gpx_buf(xc->mol->AOs.size());
	std::vector<double> gpy_buf(xc->mol->AOs.size());
	std::vector<double> gpz_buf(xc->mol->AOs.size());
	std::vector<double> temp_grad(3);
	zero_xc_data(xc);
	for(int gpt = 0; gpt < xc->g->num_gridpoints; gpt++){
		eval_bfs_grad_per_gpt(xc, phi_buf, gpx_buf, gpy_buf, gpz_buf, temp_grad, gpt);
		eval_density_grad_per_gpt(xc, phi_buf, gpx_buf, gpy_buf, gpz_buf);
		GGA_per_gpt(xc, func, phi_buf, gpx_buf, gpy_buf, gpz_buf, gpt);
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

void eval_bfs_per_gpt(XC* xc, std::vector<double>& phi_buf, int gpix){
	for(int j = 0; j < xc->mol->AOs.size(); j++){
		phi_buf[j] = xc->mol->AOs[j].evaluate(xc->g->x[gpix], xc->g->y[gpix], xc->g->z[gpix]);
	}	
}

void eval_bfs_grad_per_gpt(XC* xc, std::vector<double>& phi_buf, std::vector<double>& gpx_buf, std::vector<double>& gpy_buf, 
	std::vector<double>& gpz_buf, std::vector<double>& temp_grad, int gpix)
{
	for(int j = 0; j < xc->mol->AOs.size(); j++){
		phi_buf[j] = xc->mol->AOs[j].evaluate(xc->g->x[gpix], xc->g->y[gpix], xc->g->z[gpix]);
		temp_grad  = xc->mol->AOs[j].evaluate_gradient(xc->g->x[gpix], xc->g->y[gpix], xc->g->z[gpix]);
		gpx_buf[j] = temp_grad[0];
		gpy_buf[j] = temp_grad[1];
		gpz_buf[j] = temp_grad[2];
	}	
}

void eval_density_per_gpt(XC* xc, const std::vector<double>& phi_buf){
	if(xc->restricted){xc->rho = density(phi_buf, *(xc->P));}
	else{
		xc->rho_a = density(phi_buf, *(xc->P_A));
		xc->rho_b = density(phi_buf, *(xc->P_B));
		xc->rho = xc->rho_a + xc->rho_b;
	}
}

void eval_density_grad_per_gpt(XC* xc, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf)
{
	if(xc->restricted){
		xc->rho = density(phi_buf, *(xc->P));
		xc->gradient_rho = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc->P));
	}
	else{
		xc->rho_a = density(phi_buf, *(xc->P_A));
		xc->rho_b = density(phi_buf, *(xc->P_B));
		xc->rho = xc->rho_a + xc->rho_b;
		xc->gradient_rho_a = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc->P_A));
		xc->gradient_rho_b = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc->P_B));
		xc->gradient_rho[0] = xc->gradient_rho_a[0] + xc->gradient_rho_b[0]; 
		xc->gradient_rho[1] = xc->gradient_rho_a[1] + xc->gradient_rho_b[1]; 
		xc->gradient_rho[2] = xc->gradient_rho_a[2] + xc->gradient_rho_b[2]; 
	}
}

void LDA_per_gpt(XC* xc, LDA_ret (*func)(XC*), const std::vector<double>& phi_buf, int gpix){
	if(xc->rho < 1e-20){return;}
	std::vector<Matrix*> fxc(2);
	fxc[0] = (xc->restricted ? xc->FXC : xc->FXC_A);
	fxc[1] = (xc->restricted ? nullptr : xc->FXC_B);
	const int dim   = (xc->restricted ? xc->FXC->rows : xc->FXC_A->rows);
	const int spins = (xc->restricted ? 1 : 2);

	LDA_ret ret = func(xc);
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

void GGA_per_gpt(XC* xc, GGA_ret (*func)(XC*), const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int gpix)
{
	if(xc->rho < 1e-20){return;}
	std::vector<Matrix*> fxc(2);
	fxc[0] = (xc->restricted ? xc->FXC : xc->FXC_A);
	fxc[1] = (xc->restricted ? nullptr : xc->FXC_B);
	std::vector<std::vector<double>*> grho(2);
	grho[0] = (xc->restricted ? &xc->gradient_rho : &xc->gradient_rho_a);
	grho[1] = (xc->restricted ? nullptr           : &xc->gradient_rho_b);
	const int dim   = (xc->restricted ? xc->FXC->rows : xc->FXC_A->rows);
	const int spins = (xc->restricted ? 1 : 2);

	GGA_ret ret = func(xc);
	xc->E_XC += xc->g->w[gpix] * ret.e_XC;
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			fxc[s]->matrix[mu][mu] += xc->g->w[gpix] * (
				phi_buf[mu] * ret.drho_XC[s] * phi_buf[mu] + 
				GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, mu, s)
			);
			for(int nu = 0; nu < mu; nu++){
				fxc[s]->matrix[mu][nu] += xc->g->w[gpix] * (
					phi_buf[mu] * ret.drho_XC[s] * phi_buf[nu] + 
					GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, nu, s)
				);
				fxc[s]->matrix[nu][mu] = fxc[s]->matrix[mu][nu];
			}
		}
	}
}

// 2\sum_{\sigma'} \frac{de_{XC}^{GGA}}{d\gamma_{\sigma\sigma'}} \nabla\rho_{\sigma'} \cdot \nabla(\phi_{\mu} \phi_{\nu})
double GGA_F_second_term(const GGA_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, const std::vector<std::vector<double>*>& grho, 
	int mu, int nu, int s)
{
	double result = ret.dgamma_XC[s] * (
		phi_buf[nu] * ((*grho[s])[0] * gpx_buf[mu] + (*grho[s])[1] * gpy_buf[mu] + (*grho[s])[2] * gpz_buf[mu]) + 
		phi_buf[mu] * ((*grho[s])[0] * gpx_buf[nu] + (*grho[s])[1] * gpy_buf[nu] + (*grho[s])[2] * gpz_buf[nu])
	);
	if(ret.dgamma_XC.size()==3){
		const int s_opp = (s==0 ? 1 : 0);
		result += ret.dgamma_XC[2] * (
			phi_buf[nu] * ((*grho[s_opp])[0] * gpx_buf[mu] + (*grho[s_opp])[1] * gpy_buf[mu] + (*grho[s_opp])[2] * gpz_buf[mu]) + 
			phi_buf[mu] * ((*grho[s_opp])[0] * gpx_buf[nu] + (*grho[s_opp])[1] * gpy_buf[nu] + (*grho[s_opp])[2] * gpz_buf[nu])
		);
	}
	return 2 * result; 
}
