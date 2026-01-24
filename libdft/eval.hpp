#ifndef EVALHEADERDEF
#define EVALHEADERDEF

#include <vector>

class  XC;
struct XC_ret;

// Evaluate E_XC, F_XC
void LDA(XC* xc, XC_ret (*func)(XC*));

// Helpers
void zero_xc_data(XC* inp);
void eval_basis_funcs_per_gpt(XC* xc, std::vector<double>& phi_buf, int gpix);
void eval_density_per_gpt(XC* xc, const std::vector<double>& phi_buf);

void LDA_per_gpt(XC* xc, const std::vector<double>& phi_buf, XC_ret(*func)(XC*), int gpix);

#endif
