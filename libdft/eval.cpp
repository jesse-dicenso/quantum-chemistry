#include "eval.hpp"
#include "func.hpp"

// Evaluate E_XC, F_XC //////////////////////////////////////////

void LDA(XC* xc, LDA_ret (*func)(XC*)){
	const int size_g = xc->g->num_gridpoints;
	const int size_phi = xc->mol->AOs.size();
	zero_xc_data(xc);
	#pragma omp parallel
	{
		std::vector<double> phi_buf(size_phi);
		#pragma omp for
		for(int g = 0; g < size_g; g++){
			eval_bfs_per_gpt(xc, phi_buf, g);
			eval_density_per_gpt(xc, phi_buf);
			LDA_per_gpt(xc, func, phi_buf, g);
		}
	}
}

void GGA(XC* xc, GGA_ret (*func)(XC*)){
	const int size_g = xc->g->num_gridpoints;
	zero_xc_data(xc);
	std::vector<double> phi_buf(xc->mol->AOs.size());
	std::vector<double> gpx_buf(xc->mol->AOs.size());
	std::vector<double> gpy_buf(xc->mol->AOs.size());
	std::vector<double> gpz_buf(xc->mol->AOs.size());
	std::vector<double> temp_grad(3);
	for(int gpt = 0; gpt < size_g; gpt++){
		eval_bfs_grad_per_gpt(xc, phi_buf, gpx_buf, gpy_buf, gpz_buf, temp_grad, gpt);
		eval_density_grad_per_gpt(xc, phi_buf, gpx_buf, gpy_buf, gpz_buf);
		GGA_per_gpt(xc, func, phi_buf, gpx_buf, gpy_buf, gpz_buf, gpt);
	}
}

void MGGA(XC* xc, MGGA_ret (*func)(XC*)){
	const int size_g = xc->g->num_gridpoints;
	std::vector<double> phi_buf(xc->mol->AOs.size());
	std::vector<double> gpx_buf(xc->mol->AOs.size());
	std::vector<double> gpy_buf(xc->mol->AOs.size());
	std::vector<double> gpz_buf(xc->mol->AOs.size());
	std::vector<double> temp_grad(3);
	zero_xc_data(xc);
	int &gpt = xc->main_iter;
	for(gpt = 0; gpt < size_g; gpt++){
		eval_bfs_grad_per_gpt(xc, phi_buf, gpx_buf, gpy_buf, gpz_buf, temp_grad, gpt);
		eval_density_grad_ke_per_gpt(xc, phi_buf, gpx_buf, gpy_buf, gpz_buf);
		MGGA_per_gpt(xc, func, phi_buf, gpx_buf, gpy_buf, gpz_buf, gpt);
	}
}

// Helpers //////////////////////////////////////////////////////

void zero_xc_data(XC* xc){
	xc->main_iter = 0;
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
	const int size_p = xc->mol->AOs.size();
	const std::vector<GF>& bfs = xc->mol->AOs;
	const std::vector<double>& gx = xc->g->x;
	const std::vector<double>& gy = xc->g->y;
	const std::vector<double>& gz = xc->g->z;
	for(int j = 0; j < size_p; j++){
		phi_buf[j] = bfs[j].evaluate(gx[gpix], gy[gpix], gz[gpix]);
	}	
}

void eval_bfs_per_gpt(XC* xc, Matrix& phi_buf, int gpix){
	assert(phi_buf.cols==1);
	const int size_p = xc->mol->AOs.size();
	const std::vector<GF>& bfs = xc->mol->AOs;
	const std::vector<double>& gx = xc->g->x;
	const std::vector<double>& gy = xc->g->y;
	const std::vector<double>& gz = xc->g->z;
	for(int j = 0; j < size_p; j++){
		phi_buf.matrix[j][0] = bfs[j].evaluate(gx[gpix], gy[gpix], gz[gpix]);
	}	
}

void eval_bfs_grad_per_gpt(XC* xc, std::vector<double>& phi_buf, std::vector<double>& gpx_buf, std::vector<double>& gpy_buf, 
	std::vector<double>& gpz_buf, std::vector<double>& temp_grad, int gpix)
{
	const int size_p = xc->mol->AOs.size();
	const std::vector<GF>& bfs = xc->mol->AOs;
	const std::vector<double>& gx = xc->g->x;
	const std::vector<double>& gy = xc->g->y;
	const std::vector<double>& gz = xc->g->z;
	for(int j = 0; j < size_p; j++){
		phi_buf[j] = bfs[j].evaluate(gx[gpix], gy[gpix], gz[gpix]);
		temp_grad  = bfs[j].evaluate_gradient(gx[gpix], gy[gpix], gz[gpix]);
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

void eval_density_grad_ke_per_gpt(XC* xc, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf)
{
	if(xc->restricted){
		xc->rho = density(phi_buf, *(xc->P));
		xc->gradient_rho = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc->P));
		xc->ke_density   = ke_density(gpx_buf, gpy_buf, gpz_buf, *(xc->P));
	}
	else{
		Matrix* pa = xc->P_A;
		Matrix* pb = xc->P_B;
		xc->rho_a = density(phi_buf, *pa);
		xc->rho_b = density(phi_buf, *pb);
		xc->rho = xc->rho_a + xc->rho_b;
		xc->gradient_rho_a = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
		xc->gradient_rho_b = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
		xc->gradient_rho[0] = xc->gradient_rho_a[0] + xc->gradient_rho_b[0]; 
		xc->gradient_rho[1] = xc->gradient_rho_a[1] + xc->gradient_rho_b[1]; 
		xc->gradient_rho[2] = xc->gradient_rho_a[2] + xc->gradient_rho_b[2];
		xc->ke_density_a = ke_density(gpx_buf, gpy_buf, gpz_buf, *pa);
		xc->ke_density_b = ke_density(gpx_buf, gpy_buf, gpz_buf, *pb);
	}
}

void LDA_per_gpt(XC* xc, LDA_ret (*func)(XC*), const std::vector<double>& phi_buf, int gpix){
	if(xc->rho < 1e-20){return;}
	Matrix* fxc[2];
	fxc[0] = (xc->restricted ? xc->FXC : xc->FXC_A);
	fxc[1] = (xc->restricted ? nullptr : xc->FXC_B);
	const int dim   = (xc->restricted ? xc->FXC->rows : xc->FXC_A->rows);
	const int spins = (xc->restricted ? 1 : 2);
	const double w = xc->g->w[gpix];
	LDA_ret ret = func(xc);
	xc->E_XC += w * ret.e_XC;
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			fxc[s]->matrix[mu][mu] += w * phi_buf[mu] * ret.v_XC[s] * phi_buf[mu];
			for(int nu = 0; nu < mu; nu++){
				fxc[s]->matrix[mu][nu] += w * phi_buf[mu] * ret.v_XC[s] * phi_buf[nu];
				fxc[s]->matrix[nu][mu] = fxc[s]->matrix[mu][nu];
			}
		}
	}
}

void GGA_per_gpt(XC* xc, GGA_ret (*func)(XC*), const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int gpix)
{
	if(xc->rho < 1e-20){return;}
	Matrix* fxc[2];
	fxc[0] = (xc->restricted ? xc->FXC : xc->FXC_A);
	fxc[1] = (xc->restricted ? nullptr : xc->FXC_B);
	std::vector<double>* grho[2];
	grho[0] = (xc->restricted ? &xc->gradient_rho : &xc->gradient_rho_a);
	grho[1] = (xc->restricted ? nullptr           : &xc->gradient_rho_b);
	const int dim   = (xc->restricted ? xc->FXC->rows : xc->FXC_A->rows);
	const int spins = (xc->restricted ? 1 : 2);
	const double w = xc->g->w[gpix];
	GGA_ret ret = func(xc);
	xc->E_XC += w * ret.e_XC;
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			fxc[s]->matrix[mu][mu] += w * (
				phi_buf[mu] * ret.drho_XC[s] * phi_buf[mu] + 
				GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, mu, s)
			);
			for(int nu = 0; nu < mu; nu++){
				fxc[s]->matrix[mu][nu] += w * (
					phi_buf[mu] * ret.drho_XC[s] * phi_buf[nu] + 
					GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, nu, s)
				);
				fxc[s]->matrix[nu][mu] = fxc[s]->matrix[mu][nu];
			}
		}
	}
}

double GGA_F_second_term(const GGA_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, std::vector<double>* grho[2], int mu, int nu, int s)
{
	double result = 2 * ret.dgamma_XC[s] * (
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
	return result; 
}

void MGGA_per_gpt(XC* xc, MGGA_ret (*func)(XC*), const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int gpix)
{
	if(xc->rho < 1e-20){return;}
	Matrix* fxc[2];
	fxc[0] = (xc->restricted ? xc->FXC : xc->FXC_A);
	fxc[1] = (xc->restricted ? nullptr : xc->FXC_B);
	std::vector<double>* grho[2];
	grho[0] = (xc->restricted ? &xc->gradient_rho : &xc->gradient_rho_a);
	grho[1] = (xc->restricted ? nullptr           : &xc->gradient_rho_b);
	const int dim   = (xc->restricted ? xc->FXC->rows : xc->FXC_A->rows);
	const int spins = (xc->restricted ? 1 : 2);
	const double w = xc->g->w[gpix];
	MGGA_ret mgga_ret = func(xc);
	GGA_ret tmp_gga_ret;
	tmp_gga_ret.drho_XC = mgga_ret.drho_XC;
	tmp_gga_ret.dgamma_XC = mgga_ret.dgamma_XC;
	xc->E_XC += w * mgga_ret.e_XC;
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			fxc[s]->matrix[mu][mu] += w * (
				phi_buf[mu] * mgga_ret.drho_XC[s] * phi_buf[mu] + 
				GGA_F_second_term(tmp_gga_ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, mu, s) + 
				MGGA_F_third_term(mgga_ret, gpx_buf, gpy_buf, gpz_buf, mu, mu, s)
			);
			for(int nu = 0; nu < mu; nu++){
				fxc[s]->matrix[mu][nu] += w * (
					phi_buf[mu] * mgga_ret.drho_XC[s] * phi_buf[nu] + 
					GGA_F_second_term(tmp_gga_ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, nu, s) +
					MGGA_F_third_term(mgga_ret, gpx_buf, gpy_buf, gpz_buf, mu, nu, s)
				);
				fxc[s]->matrix[nu][mu] = fxc[s]->matrix[mu][nu];
			}
		}
	}
}

// FDO potential: \frac{1}{2} \frac{de^{MGGA}_{XC}}{d\tau_{\sigma}} \nabla \phi_{\mu} \cdot \nabla \phi_{\nu}
double MGGA_F_third_term(const MGGA_ret& ret, const std::vector<double>& gpx_buf, const std::vector<double>& gpy_buf, 
	const std::vector<double>& gpz_buf, int mu, int nu, int s)
{
	return ret.dtau_XC[s] * (gpx_buf[mu] * gpx_buf[nu] + gpy_buf[mu] * gpy_buf[nu] + gpz_buf[mu] * gpz_buf[nu]);
}

