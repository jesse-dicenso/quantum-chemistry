#include "eval.hpp"
#include "func.hpp"

// Evaluate E_XC, F_XC //////////////////////////////////////////

void scf_xc_call(XC* xc){
    if(xc->isHF){HFX(xc);}
    else if(xc->isHFSN){HFSNX(xc);}
    else if(xc->isLDA){LDA(xc);}
    else if((xc->isGGA) || (xc->isMGGA)){GGA_MGGA(xc);}
    if(xc->isNLC){xc->nlc_functional(xc);}
}

void LDA(XC* xc){
	const int spins = (xc->restricted ? 1 : 2);
    const int mat_dim = xc->F_XC[0]->rows;
    double E_XC = 0.0;
	zero_xc_data(xc, spins);
	#pragma omp parallel
	{
        XC_inp inp;
        XC_ret ret;
        mat_alloc(ret.F_XC, spins, mat_dim, mat_dim);
		std::vector<double> phi_buf(xc->mol->AOs.size());
		#pragma omp for reduction(+:E_XC)
		for(int g = 0; g < xc->g->num_gridpoints; g++){
			eval_bfs_per_gpt(*xc, phi_buf, g);
			eval_density_per_gpt(*xc, inp, phi_buf);
			E_XC += LDA_per_gpt(*xc, inp, ret, phi_buf, spins, g);
		}
        #pragma omp critical
        { 
            for(int s = 0; s < spins; s++){
                *(xc->F_XC[s]) = *(xc->F_XC[s]) + ret.F_XC[s];
            }
        }
	}
    xc->E_XC = E_XC;
}

void GGA_MGGA(XC* xc){
	const int spins = (xc->restricted ? 1 : 2);
    const int mat_dim = xc->F_XC[0]->rows;
    double E_XC = 0.0;
	zero_xc_data(xc, spins);

    double (*dft_per_gpt)(const XC&, const XC_inp&, XC_ret&, const std::vector<double>&, const std::vector<double>&, 
        const std::vector<double>&, const std::vector<double>&, int, int);
    void (*eval_properties)(const XC&, XC_inp&, const std::vector<double>&, const std::vector<double>&, 
            const std::vector<double>&, const std::vector<double>&);
    if(xc->isGGA){
        dft_per_gpt = GGA_per_gpt;
        eval_properties = eval_density_grad_per_gpt;
    }
    else if(xc->isMGGA){
        dft_per_gpt = MGGA_per_gpt;
        eval_properties = eval_density_grad_ke_per_gpt;
    }
    else{assert((xc->isGGA) || (xc->isMGGA));}
	#pragma omp parallel
	{
        XC_inp inp;
        XC_ret ret;
        mat_alloc(ret.F_XC, spins, mat_dim, mat_dim);
		std::vector<double> phi_buf(xc->mol->AOs.size());
        std::vector<double> gpx_buf(xc->mol->AOs.size());
        std::vector<double> gpy_buf(xc->mol->AOs.size());
        std::vector<double> gpz_buf(xc->mol->AOs.size());
        std::vector<double> temp_grad(3);
		#pragma omp for reduction(+:E_XC)
		for(int g = 0; g < xc->g->num_gridpoints; g++){
		    eval_bfs_grad_per_gpt(*xc, phi_buf, gpx_buf, gpy_buf, gpz_buf, temp_grad, g);
		    eval_properties(*xc, inp, phi_buf, gpx_buf, gpy_buf, gpz_buf);
			E_XC += dft_per_gpt(*xc, inp, ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, spins, g);
		}
        #pragma omp critical
        { 
            for(int s = 0; s < spins; s++){
                *(xc->F_XC[s]) = *(xc->F_XC[s]) + ret.F_XC[s];
            }
        }
	}
    xc->E_XC = E_XC;
}

// Helpers //////////////////////////////////////////////////////

void zero_xc_data(XC* xc, int spins){
	xc->E_XC = 0.0;
    for(int s = 0; s < spins; s++){
	    assert(xc->F_XC[s]!=nullptr);
		*(xc->F_XC[s]) = zero(xc->F_XC[s]->rows, xc->F_XC[s]->cols);
    }
}

void eval_bfs_per_gpt(const XC& xc, std::vector<double>& phi_buf, int gpix){
	const int size_p = xc.mol->AOs.size();
	const std::vector<GF>& bfs = xc.mol->AOs;
	const std::vector<double>& gx = xc.g->x;
	const std::vector<double>& gy = xc.g->y;
	const std::vector<double>& gz = xc.g->z;
	for(int j = 0; j < size_p; j++){
		phi_buf[j] = bfs[j].evaluate(gx[gpix], gy[gpix], gz[gpix]);
	}	
}

void eval_bfs_per_gpt(const XC& xc, Matrix& phi_buf, int gpix){
	assert(phi_buf.cols==1);
	const int size_p = xc.mol->AOs.size();
	const std::vector<GF>& bfs = xc.mol->AOs;
	const std::vector<double>& gx = xc.g->x;
	const std::vector<double>& gy = xc.g->y;
	const std::vector<double>& gz = xc.g->z;
	for(int j = 0; j < size_p; j++){
		phi_buf(j, 0) = bfs[j].evaluate(gx[gpix], gy[gpix], gz[gpix]);
	}	
}

void eval_bfs_per_batch(const XC& xc, Matrix& phi_buf, int g_start, int g_end){
	assert(phi_buf.cols==(g_start - g_end));
	const int size_p = xc.mol->AOs.size();
	const std::vector<GF>& bfs = xc.mol->AOs;
	const std::vector<double>& gx = xc.g->x;
	const std::vector<double>& gy = xc.g->y;
	const std::vector<double>& gz = xc.g->z;
	for(int j = 0; j < size_p; j++){
        for(int sub_g = g_start; sub_g < g_end; sub_g++){
		    phi_buf(j, sub_g) = bfs[j].evaluate(gx[sub_g], gy[sub_g], gz[sub_g]);
        }
	}	
}

void eval_bfs_grad_per_gpt(const XC& xc, std::vector<double>& phi_buf, std::vector<double>& gpx_buf, std::vector<double>& gpy_buf, 
	std::vector<double>& gpz_buf, std::vector<double>& temp_grad, int gpix)
{
	const int size_p = xc.mol->AOs.size();
	const std::vector<GF>& bfs = xc.mol->AOs;
	const std::vector<double>& gx = xc.g->x;
	const std::vector<double>& gy = xc.g->y;
	const std::vector<double>& gz = xc.g->z;
	for(int j = 0; j < size_p; j++){
		phi_buf[j] = bfs[j].evaluate(gx[gpix], gy[gpix], gz[gpix]);
		temp_grad  = bfs[j].evaluate_gradient(gx[gpix], gy[gpix], gz[gpix]);
		gpx_buf[j] = temp_grad[0];
		gpy_buf[j] = temp_grad[1];
		gpz_buf[j] = temp_grad[2];
	}	
}

void eval_density_per_gpt(const XC& xc, XC_inp& inp, const std::vector<double>& phi_buf){
	if(xc.restricted){inp.rho = density(phi_buf, *(xc.P[0]));}
	else{
		inp.rho_a = density(phi_buf, *(xc.P[0]));
		inp.rho_b = density(phi_buf, *(xc.P[1]));
		inp.rho = inp.rho_a + inp.rho_b;
	}
}

void eval_density_grad_per_gpt(const XC& xc, XC_inp& inp, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf)
{
	if(xc.restricted){
		inp.rho = density(phi_buf, *(xc.P[0]));
		inp.gradient_rho = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc.P[0]));
	}
	else{
		inp.rho_a = density(phi_buf, *(xc.P[0]));
		inp.rho_b = density(phi_buf, *(xc.P[1]));
		inp.rho = inp.rho_a + inp.rho_b;
		inp.gradient_rho_a = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc.P[0]));
		inp.gradient_rho_b = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc.P[1]));
		inp.gradient_rho[0] = inp.gradient_rho_a[0] + inp.gradient_rho_b[0]; 
		inp.gradient_rho[1] = inp.gradient_rho_a[1] + inp.gradient_rho_b[1]; 
		inp.gradient_rho[2] = inp.gradient_rho_a[2] + inp.gradient_rho_b[2]; 
	}
}

void eval_density_grad_ke_per_gpt(const XC& xc, XC_inp& inp, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf)
{
	if(xc.restricted){
		inp.rho = density(phi_buf, *(xc.P[0]));
		inp.gradient_rho = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *(xc.P[0]));
		inp.ke_density   = ke_density(gpx_buf, gpy_buf, gpz_buf, *(xc.P[0]));
	}
	else{
		Matrix* pa = xc.P[0];
		Matrix* pb = xc.P[1];
		inp.rho_a = density(phi_buf, *pa);
		inp.rho_b = density(phi_buf, *pb);
		inp.rho = inp.rho_a + inp.rho_b;
		inp.gradient_rho_a = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
		inp.gradient_rho_b = density_gradient(phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
		inp.gradient_rho[0] = inp.gradient_rho_a[0] + inp.gradient_rho_b[0]; 
		inp.gradient_rho[1] = inp.gradient_rho_a[1] + inp.gradient_rho_b[1]; 
		inp.gradient_rho[2] = inp.gradient_rho_a[2] + inp.gradient_rho_b[2];
		inp.ke_density_a = ke_density(gpx_buf, gpy_buf, gpz_buf, *pa);
		inp.ke_density_b = ke_density(gpx_buf, gpy_buf, gpz_buf, *pb);
	}
}

// return e_XC * w[gpix]
double LDA_per_gpt(const XC& xc, const XC_inp& inp, XC_ret& ret, const std::vector<double>& phi_buf, int spins, int gpix){
	if(inp.rho < 1e-20){return 0.0;}
	const int dim = ret.F_XC[0].rows;
    const double w = xc.g->w[gpix];
    std::vector<Matrix>& F_XC = ret.F_XC;
	xc.xc_functional(xc, inp, ret);
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			F_XC[s](mu, mu) += w * phi_buf[mu] * ret.drho_XC[s] * phi_buf[mu];
			for(int nu = 0; nu < mu; nu++){
				F_XC[s](mu, nu) += w * phi_buf[mu] * ret.drho_XC[s] * phi_buf[nu];
				F_XC[s](nu, mu) = F_XC[s](mu, nu);
			}
		}
	}
    return w * ret.e_XC;
}

double GGA_per_gpt(const XC& xc, const XC_inp& inp, XC_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int spins, int gpix)
{
	if(inp.rho < 1e-20){return 0.0;}
	const std::vector<double>* grho[2] = {
        (xc.restricted ? &inp.gradient_rho : &inp.gradient_rho_a),
    	(xc.restricted ? nullptr           : &inp.gradient_rho_b)
    };
	const int dim = ret.F_XC[0].rows;
	const double w = xc.g->w[gpix];
    std::vector<Matrix>& F_XC = ret.F_XC;
	xc.xc_functional(xc, inp, ret);
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			F_XC[s](mu, mu) += w * (
				phi_buf[mu] * ret.drho_XC[s] * phi_buf[mu] + 
				GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, mu, s)
			);
			for(int nu = 0; nu < mu; nu++){
				F_XC[s](mu, nu) += w * (
					phi_buf[mu] * ret.drho_XC[s] * phi_buf[nu] + 
					GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, nu, s)
				);
				F_XC[s](nu, mu) = F_XC[s](mu, nu);
			}
		}
	}
    return w * ret.e_XC;
}

double MGGA_per_gpt(const XC& xc, const XC_inp& inp, XC_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, int spins, int gpix)
{
	if(inp.rho < 1e-20){return 0.0;}
	const std::vector<double>* grho[2] = {
        (xc.restricted ? &inp.gradient_rho : &inp.gradient_rho_a),
    	(xc.restricted ? nullptr           : &inp.gradient_rho_b)
    };
	const int dim = ret.F_XC[0].rows;
	const double w = xc.g->w[gpix];
    std::vector<Matrix>& F_XC = ret.F_XC;
	xc.xc_functional(xc, inp, ret);
	for(int s = 0; s < spins; s++){
		for(int mu = 0; mu < dim; mu++){
			F_XC[s](mu, mu) += w * (
				phi_buf[mu] * ret.drho_XC[s] * phi_buf[mu] + 
				GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, mu, s) +
				MGGA_F_third_term(ret, gpx_buf, gpy_buf, gpz_buf, mu, mu, s)
			);
			for(int nu = 0; nu < mu; nu++){
				F_XC[s](mu, nu) += w * (
					phi_buf[mu] * ret.drho_XC[s] * phi_buf[nu] + 
					GGA_F_second_term(ret, phi_buf, gpx_buf, gpy_buf, gpz_buf, grho, mu, nu, s) +
				    MGGA_F_third_term(ret, gpx_buf, gpy_buf, gpz_buf, mu, nu, s)
				);
				F_XC[s](nu, mu) = F_XC[s](mu, nu);
			}
		}
	}
    return w * ret.e_XC;
}

double GGA_F_second_term(const XC_ret& ret, const std::vector<double>& phi_buf, const std::vector<double>& gpx_buf, 
	const std::vector<double>& gpy_buf, const std::vector<double>& gpz_buf, const std::vector<double>* grho[2], int mu, int nu, int s)
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

// FDO-type MGGA Fock matrix term
double MGGA_F_third_term(const XC_ret& ret, const std::vector<double>& gpx_buf, const std::vector<double>& gpy_buf, 
	const std::vector<double>& gpz_buf, int mu, int nu, int s)
{
	return ret.dtau_XC[s] * (gpx_buf[mu] * gpx_buf[nu] + gpy_buf[mu] * gpy_buf[nu] + gpz_buf[mu] * gpz_buf[nu]);
}

