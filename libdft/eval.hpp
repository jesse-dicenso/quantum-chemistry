#ifndef EVALHEADERDEF
#define EVALHEADERDEF

#include<vector>

class  XC;
struct XC_ret;

void zero_xc_data(XC* inp);
void eval_basis_funcs_per_gpt(XC* xc, std::vector<double>& phi_buf, int gpix);
void eval_density_per_gpt(XC* xc, const std::vector<double>& phi_buf);

void LDA_per_gpt(XC* xc, const std::vector<double>& phi_buf, XC_ret(XC*) func, int gpix);

// LDA //
void LDA(XC* xc, XC_ret(XC*) func){
	int num_gpts = xc->g->num_gridpoints;
	std::vector<double> phi_buf(xc->mol->AOs.size());
	zero_xc_data(xc);
	for(int i = 0; i < num_gpts; i++){
		eval_basis_funcs_per_gpt(xc, phi_buf, i);
		eval_density_per_gpt(xc, phi_buf);
		LDA_per_gpt(xc, phi_buf, func);
	}
}

#endif
