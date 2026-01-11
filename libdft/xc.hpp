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
template <int sflag, typename F, typename... Args>
XC_ret F_XC_LDA(const XC_inp& inp, F&& v_LDA, const Args&... args){
	Molecule* m = inp.mol;
	grid*     g = inp.g;
	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(m->AOs.size());
	// Restricted: sflag = 0
	if constexpr (sflag==0){
		Matrix* p = inp.PT;
		Matrix F_XC(p->rows, p->cols);
		Matrix null;
		double rho;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho = density(g->x[i], g->y[i], g->z[i], phi_buf, *p);
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
	else if constexpr (sflag==1){
		Matrix* pa = inp.PA;
		Matrix* pb = inp.PB;
		Matrix F_XC_A(pa->rows, pa->cols);
		Matrix F_XC_B(pb->rows, pb->cols);
		double rho_a, rho_b;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho_a = density(g->x[i], g->y[i], g->z[i], phi_buf, *pa);
			rho_b = density(g->x[i], g->y[i], g->z[i], phi_buf, *pb);
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
	else if constexpr(sflag==2){
		Matrix* pa = inp.PA;
		Matrix* pb = inp.PB;
		Matrix F_XC_A(pa->rows, pa->cols);
		Matrix F_XC_B(pb->rows, pb->cols);
		double rho_a, rho_b;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho_a = density(g->x[i], g->y[i], g->z[i], phi_buf, *pa);
			rho_b = density(g->x[i], g->y[i], g->z[i], phi_buf, *pb);
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

template <int sflag, typename F, typename... Args>
double E_XC_LDA(const XC_inp& inp, F&& e_LDA, const Args&... args){
	Molecule* m = inp.mol;
	grid*     g = inp.g;
	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(m->AOs.size());
	double E_XC = 0;
	// Restricted: sflag = 0
	if constexpr (sflag==0){
		Matrix* p = inp.PT;
		double rho;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho = density(g->x[i], g->y[i], g->z[i], phi_buf, *p);
			E_XC += g->w[i] * e_LDA(rho, args...);
		}
	}
	// Unrestricted: sflag = 1
	else if constexpr (sflag==1){
		Matrix* pa = inp.PA;
		Matrix* pb = inp.PB;
		double rho_a, rho_b;
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
			}
			rho_a = density(g->x[i], g->y[i], g->z[i], phi_buf, *pa);
			rho_b = density(g->x[i], g->y[i], g->z[i], phi_buf, *pb);
			E_XC += g->w[i] * e_LDA(rho_a, rho_b, args...);
		}
	}
	else{assert((sflag==0) || (sflag==1));}
	return E_XC;
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

XC_ret R_PW92_c(const XC_inp& inp);
double R_PW92_c_E(const XC_inp& inp);
XC_ret U_PW92_c(const XC_inp& inp);
double U_PW92_c_E(const XC_inp& inp);
XC_ret R_PW92(const XC_inp& inp);
double R_PW92_E(const XC_inp& inp);
XC_ret U_PW92(const XC_inp& inp);
double U_PW92_E(const XC_inp& inp);

// GGA //
template <int sflag, typename F>
XC_ret F_XC_GGA(const XC_inp& inp, F&& v_GGA){
	Molecule* m = inp.mol;
	grid*     g = inp.g;
	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(m->AOs.size());
	std::vector<double> gpx_buf(m->AOs.size());
	std::vector<double> gpy_buf(m->AOs.size());
	std::vector<double> gpz_buf(m->AOs.size());
	std::vector<double> tmp_grd(3);
	// Restricted: sflag = 0
	if constexpr (sflag==0){
		Matrix* p = inp.PT;
		Matrix F_XC(p->rows, p->cols);
		Matrix null;
		double rho;
		std::vector<double> grho(3);
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
				tmp_grd = m->AOs[j].evaluate_gradient(g->x[i], g->y[i], g->z[i]);
				gpx_buf[j] = tmp_grd[0];
				gpy_buf[j] = tmp_grd[1];
				gpz_buf[j] = tmp_grd[2];
			}
			rho  = density(g->x[i], g->y[i], g->z[i], phi_buf, *p);
			grho = density_gradient(g->x[i], g->y[i], g->z[i], phi_buf, gpx_buf, gpy_buf, gpz_buf, *p);
			for(int mu = 0; mu < F_XC.rows; mu++){
				F_XC.matrix[mu][mu] += g->w[i] * v_GGA(rho, grho, phi_buf[mu], phi_buf[mu], 
													   gpx_buf[mu], gpy_buf[mu], gpz_buf[mu],
													   gpx_buf[mu], gpy_buf[mu], gpz_buf[mu]);
				for(int nu = 0; nu < mu; nu++){
					F_XC.matrix[mu][mu] += g->w[i] * v_GGA(rho, grho, phi_buf[mu], phi_buf[nu], 
														   gpx_buf[mu], gpy_buf[mu], gpz_buf[mu],
														   gpx_buf[nu], gpy_buf[nu], gpz_buf[nu]);
					F_XC.matrix[nu][mu] = F_XC.matrix[mu][nu];
				}
			}
		}
		return {F_XC, null};
	}
	// Unrestricted: sflag = 1
	else if constexpr (sflag==1){
		Matrix* pa = inp.PA;
		Matrix* pb = inp.PB;
		Matrix F_XC_A(pa->rows, pa->cols);
		Matrix F_XC_B(pb->rows, pb->cols);
		double rho_a, rho_b;
		std::vector<double> grho_a(3);
		std::vector<double> grho_b(3);
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
				tmp_grd = m->AOs[j].evaluate_gradient(g->x[i], g->y[i], g->z[i]);
				gpx_buf[j] = tmp_grd[0];
				gpy_buf[j] = tmp_grd[1];
				gpz_buf[j] = tmp_grd[2];
			}
			rho_a = density(g->x[i], g->y[i], g->z[i], phi_buf, *pa);
			rho_b = density(g->x[i], g->y[i], g->z[i], phi_buf, *pb);
			grho_a = density_gradient(g->x[i], g->y[i], g->z[i], phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
			grho_b = density_gradient(g->x[i], g->y[i], g->z[i], phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
			for(int mu = 0; mu < F_XC_A.rows; mu++){
				F_XC_A.matrix[mu][mu] += g->w[i] * v_GGA(rho_a, rho_b, grho_a, grho_b, phi_buf[mu], phi_buf[mu], 
													     gpx_buf[mu], gpy_buf[mu], gpz_buf[mu],
													     gpx_buf[mu], gpy_buf[mu], gpz_buf[mu], 0);
				F_XC_B.matrix[mu][mu] += g->w[i] * v_GGA(rho_a, rho_b, grho_a, grho_b, phi_buf[mu], phi_buf[mu], 
													     gpx_buf[mu], gpy_buf[mu], gpz_buf[mu],
													     gpx_buf[mu], gpy_buf[mu], gpz_buf[mu], 1);
				for(int nu = 0; nu < mu; nu++){
					F_XC_A.matrix[mu][mu] += g->w[i] * v_GGA(rho_a, rho_b, grho_a, grho_b, phi_buf[mu], phi_buf[nu], 
														     gpx_buf[mu], gpy_buf[mu], gpz_buf[mu],
														     gpx_buf[nu], gpy_buf[nu], gpz_buf[nu], 0);
					F_XC_A.matrix[nu][mu] = F_XC_A.matrix[mu][nu];
					F_XC_B.matrix[mu][mu] += g->w[i] * v_GGA(rho_a, rho_b, grho_a, grho_b, phi_buf[mu], phi_buf[nu], 
														     gpx_buf[mu], gpy_buf[mu], gpz_buf[mu],
														     gpx_buf[nu], gpy_buf[nu], gpz_buf[nu], 1);
					F_XC_B.matrix[nu][mu] = F_XC_B.matrix[mu][nu];
				}
			}
		}
		return {F_XC_A, F_XC_B};
	}
}

template <int sflag, typename F>
double E_XC_GGA(const XC_inp& inp, F&& e_GGA){
	Molecule* m = inp.mol;
	grid*     g = inp.g;
	int num_gpts = g->num_gridpoints;
	std::vector<double> phi_buf(m->AOs.size());
	std::vector<double> gpx_buf(m->AOs.size());
	std::vector<double> gpy_buf(m->AOs.size());
	std::vector<double> gpz_buf(m->AOs.size());
	std::vector<double> tmp_grd(3);
	double E_XC = 0;
	// Restricted: sflag = 0
	if constexpr (sflag==0){
		Matrix* p = inp.PT;
		double rho;
		std::vector<double> grho(3);
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
				tmp_grd = m->AOs[j].evaluate_gradient(g->x[i], g->y[i], g->z[i]);
				gpx_buf[j] = tmp_grd[0];
				gpy_buf[j] = tmp_grd[1];
				gpz_buf[j] = tmp_grd[2];
			}
			rho = density(g->x[i], g->y[i], g->z[i], phi_buf, *p);
			grho = density_gradient(g->x[i], g->y[i], g->z[i], phi_buf, gpx_buf, gpy_buf, gpz_buf, *p);
			E_XC += g->w[i] * e_GGA(rho, grho);
		}
	}
	// Unrestricted: sflag = 1
	else if constexpr (sflag==1){
		Matrix* pa = inp.PA;
		Matrix* pb = inp.PB;
		double rho_a, rho_b;
		std::vector<double> grho_a(3);
		std::vector<double> grho_b(3);
		for(int i = 0; i < num_gpts; i++){
			for(int j = 0; j < m->AOs.size(); j++){
				phi_buf[j] = m->AOs[j].evaluate(g->x[i], g->y[i], g->z[i]);
				tmp_grd = m->AOs[j].evaluate_gradient(g->x[i], g->y[i], g->z[i]);
				gpx_buf[j] = tmp_grd[0];
				gpy_buf[j] = tmp_grd[1];
				gpz_buf[j] = tmp_grd[2];
			}
			rho_a = density(g->x[i], g->y[i], g->z[i], phi_buf, *pa);
			rho_b = density(g->x[i], g->y[i], g->z[i], phi_buf, *pb);
			grho_a = density_gradient(g->x[i], g->y[i], g->z[i], phi_buf, gpx_buf, gpy_buf, gpz_buf, *pa);
			grho_b = density_gradient(g->x[i], g->y[i], g->z[i], phi_buf, gpx_buf, gpy_buf, gpz_buf, *pb);
			E_XC += g->w[i] * e_GGA(rho_a, rho_b, grho_a, grho_b);
		}
	}
	else{assert((sflag==0) || (sflag==1));}
	return E_XC;
}

XC_ret R_PBE_X(const XC_inp& inp);
double R_PBE_X_E(const XC_inp& inp);
XC_ret U_PBE_X(const XC_inp& inp);
double U_PBE_X_E(const XC_inp& inp);
XC_ret R_PBE_c(const XC_inp& inp);
double R_PBE_c_E(const XC_inp& inp);
XC_ret U_PBE_c(const XC_inp& inp);
double U_PBE_c_E(const XC_inp& inp);
XC_ret R_PBE(const XC_inp& inp);
double R_PBE_E(const XC_inp& inp);
XC_ret U_PBE(const XC_inp& inp);
double U_PBE_E(const XC_inp& inp);

#endif
