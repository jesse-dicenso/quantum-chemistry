#ifndef XCHEADERDEF
#define XCHEADERDEF

#include "dft_helper.hpp"

#include <string>
#include <functional>
#include <unordered_map>


class XC_inp{
	public:
		XC_inp(const std::string& method_name);
		
		std::string method;
	
		bool is_HF;
		bool is_LDA;
		bool is_GGA;
		
		Matrix* PT = nullptr;
		Matrix* PA = nullptr;
		Matrix* PB = nullptr;

		std::vector<std::vector<std::vector<std::vector<double>>>>* eris = nullptr;
		Molecule* mol = nullptr;
		grid* g = nullptr;
};

struct XC_ret{
	Matrix F_XC_1;	// F_XC if restricted; F_XC_a if unrestricted
	Matrix F_XC_2;	// F_XC_b if unrestricted
};

extern std::unordered_map<std::string, std::function<XC_ret(const XC_inp&)>> xc_v_register;
XC_ret F_XC(XC_inp* inp);

extern std::unordered_map<std::string, std::function<double(const XC_inp&)>> xc_E_register;
double E_XC(XC_inp* inp);

// HF //
XC_ret R_HF_X(const XC_inp& inp);
XC_ret U_HF_X(const XC_inp& inp);

// LDA //
template <typename F, typename... Args>
XC_ret F_XC_LDA(const XC_inp& inp, int sflag, F&& v_LDA, const Args&... args){
	Molecule* m = inp.mol;
	grid*     g = inp.g;
	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(num_gpts);
	// Restricted: sflag = 0
	if(sflag==0){
		Matrix* p = inp.PT;
		Matrix F_XC(p->rows, p->cols);
		Matrix null;
		double rho;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho = density2(g->x[i], g->y[i], g->z[i], phi_buf, *p);
			for(int mu = 0; mu < F_XC.rows; mu++){
				F_XC.matrix[mu][mu] += g->w[i] * phi_buf[mu] * v_LDA(rho, args...) * phi_buf[mu];
				for(int nu = 0; nu < mu; nu++){
					F_XC.matrix[mu][nu] += g->w[i] * phi_buf[mu] * v_LDA(rho, args...) * phi_buf[nu];
					F_XC.matrix[nu][mu] = F_XC.matrix[mu][nu];
				}
			}
		}
		return {F_XC, null};
	}
	// Unrestricted, rho_a, rho_b not needed together: sflag = 1
	else if(sflag==1){
		Matrix* pa = inp.PA;
		Matrix* pb = inp.PB;
		Matrix F_XC_A(pa->rows, pa->cols);
		Matrix F_XC_B(pb->rows, pb->cols);
		double rho_a, rho_b;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho_a = density2(g->x[i], g->y[i], g->z[i], phi_buf, *pa);
			rho_b = density2(g->x[i], g->y[i], g->z[i], phi_buf, *pb);
			for(int mu = 0; mu < F_XC_A.rows; mu++){
				F_XC_A.matrix[mu][mu] += g->w[i] * phi_buf[mu] * v_LDA(rho_a, args...) * phi_buf[mu];
				F_XC_B.matrix[mu][mu] += g->w[i] * phi_buf[mu] * v_LDA(rho_b, args...) * phi_buf[mu];
				for(int nu = 0; nu < mu; nu++){
					F_XC_A.matrix[mu][nu] += g->w[i] * phi_buf[mu] * v_LDA(rho_a, args...) * phi_buf[nu];
					F_XC_A.matrix[nu][mu] = F_XC_A.matrix[mu][nu];
					F_XC_B.matrix[mu][nu] += g->w[i] * phi_buf[mu] * v_LDA(rho_b, args...) * phi_buf[nu];
					F_XC_B.matrix[nu][mu] = F_XC_B.matrix[mu][nu];
				}
			}
		}
		return {F_XC_A, F_XC_B};
	}
	// Unrestricted, rho_a, rho_b needed together: sflag = 2
	else if(sflag==2){
		Matrix* pa = inp.PA;
		Matrix* pb = inp.PB;
		Matrix F_XC_A(pa->rows, pa->cols);
		Matrix F_XC_B(pb->rows, pb->cols);
		double rho_a, rho_b;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho_a = density2(g->x[i], g->y[i], g->z[i], phi_buf, *pa);
			rho_b = density2(g->x[i], g->y[i], g->z[i], phi_buf, *pb);
			for(int mu = 0; mu < F_XC_A.rows; mu++){
				F_XC_A.matrix[mu][mu] += g->w[i] * phi_buf[mu] * v_LDA(rho_a, rho_b, 0, args...) * phi_buf[mu];
				F_XC_B.matrix[mu][mu] += g->w[i] * phi_buf[mu] * v_LDA(rho_a, rho_b, 1, args...) * phi_buf[mu];
				for(int nu = 0; nu < mu; nu++){
					F_XC_A.matrix[mu][nu] += g->w[i] * phi_buf[mu] * v_LDA(rho_a, rho_b, 0, args...) * phi_buf[nu];
					F_XC_A.matrix[nu][mu] = F_XC_A.matrix[mu][nu];
					F_XC_B.matrix[mu][nu] += g->w[i] * phi_buf[mu] * v_LDA(rho_a, rho_b, 1, args...) * phi_buf[nu];
					F_XC_B.matrix[nu][mu] = F_XC_B.matrix[mu][nu];
				}
			}
		}
		return {F_XC_A, F_XC_B};
	}
	else{assert((sflag==0) || (sflag==1) || (sflag==2));}
}

XC_ret R_Slater_X(const XC_inp& inp);
double R_Slater_X_E(const XC_inp& inp);
XC_ret U_Slater_X(const XC_inp& inp);
double U_Slater_X_E(const XC_inp& inp);

XC_ret R_VWN5_c(const XC_inp& inp);
double R_VWN5_c_E(const XC_inp& inp);
XC_ret U_VWN5_c(const XC_inp& inp);
double U_VWN5_c_E(const XC_inp& inp);
XC_ret R_VWN5(const XC_inp& inp);
double R_VWN5_E(const XC_inp& inp);
XC_ret U_VWN5(const XC_inp& inp);
double U_VWN5_E(const XC_inp& inp);

// GGA //

#endif
