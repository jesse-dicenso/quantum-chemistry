#include "func.hpp"
#include "eval.hpp"

XC::XC(const std::string& method){
	if		(method.substr(0,2)=="R_"){restricted=true;}
	else if (method.substr(0,2)=="U_"){restricted=false;}
	else{throw std::invalid_argument("ERR: method must be restricted 'R_' or unrestricted 'U_'");}
    const Functional func = xc_register.at(method.substr(2));
    switch(func){
        case Functional::HF : 
            isHF = true;
            break;
        case Functional::SNX : 
            isSNX = true;
            break;
        case Functional::SLATER : 
            xc_functional = Slater;
            isLDA = true;
            break;
        case Functional::VWN5 : 
            xc_functional = VWN5;
            isLDA = true;
            break;
        case Functional::PW92 : 
            xc_functional = PW92;
            isLDA = true;
            break;
        case Functional::PBE_X : 
            xc_functional = PBE_X;
            isGGA = true;
            break;
        case Functional::PBE : 
            xc_functional = PBE;
            isGGA = true;
            break;
        case Functional::B97M_V : 
            xc_functional = B97M_V;
            nlc_functional = VV10;
            nlc_params = {6.0, 0.01};
            isMGGA = true;
            isNLC = true;
            break;
        default :
            throw std::invalid_argument("ERR: method unknown");
            break;
    }
}

const std::unordered_map<std::string, Functional> xc_register = 
{
	{ "HF"     , Functional::HF     },
	{ "SNX"    , Functional::SNX    },
	{ "Slater" , Functional::SLATER },
	{ "VWN5"   , Functional::VWN5   },
	{ "PW92"   , Functional::PW92   },
	{ "PBE_X"  , Functional::PBE_X  },
	{ "PBE"    , Functional::PBE    },
	{ "B97M-V" , Functional::B97M_V }
};

///////////////////////////////////////////////////////////////////
// HF /////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

void HFX(XC* xc){
	assert(xc->eris!=nullptr);
    const int spins = (xc->restricted ? 1 : 2);
    for(int s = 0; s < spins; s++){
        assert((xc->P[s]!=nullptr) && (xc->F_XC[s]!=nullptr));
    }
    const int dim = xc->F_XC[0]->rows;
    const double spin_factor = (xc->restricted ? 0.5 : 1.0);
    std::vector<double> fxc(spins, 0.0);
	xc->E_XC = 0.0;
	for(int mu = 0; mu < dim; mu++){
		for(int nu = 0; nu < dim; nu++){
            for(int s = 0; s < spins; s++){
                fxc[s] = 0.0;
            }
			for(int ld = 0; ld < dim; ld++){
				for(int sg = 0; sg < dim; sg++){
                    for(int s = 0; s < spins; s++){
					    fxc[s] -= (*xc->P[s])(ld, sg) * (*xc->eris)[mu][ld][sg][nu];
                    }
				}
			}
            for(int s = 0; s < spins; s++){
			    (*xc->F_XC[s])(mu, nu) = spin_factor * fxc[s]; 
			    xc->E_XC += (*xc->F_XC[s])(mu, nu) * (*xc->P[s])(mu, nu);
            }
		}
	}
	xc->E_XC *= 0.5;
}

void SNX(XC* xc){
    //auto start_K = std::chrono::high_resolution_clock::now();
    assert((xc->g!=nullptr) && (xc->mol!=nullptr));
    
    const double snx_int_thresh = xc->snx_thresh_k;
    const Matrix& Vs = xc->snx_screen;
    
    const int spins = (xc->restricted ? 1 : 2);
    const double spin_factor = (xc->restricted ? 0.5 : 1.0);
    for(int s = 0; s < spins; s++){
        assert((xc->P[s]!=nullptr) && (xc->F_XC[s]!=nullptr));
    }
	const std::vector<Matrix*>& p = xc->P;
    const int size_p = xc->mol->AOs.size();
    
    const int size_g = xc->g->num_gridpoints;
    constexpr int target_batch_size = 64;
    const int num_batches = (size_g + (target_batch_size - 1)) / target_batch_size; // round up int div	
   
    zero_xc_data(xc, spins);
    #pragma omp parallel
    {
	    Matrix X;
        Tensor3 A;
        std::vector<Matrix> G_T;
        std::vector<Matrix> fxc;
        mat_alloc(fxc, spins, size_p, size_p);
        #pragma omp for schedule(dynamic)
        for(int gb = 0; gb < num_batches; gb++){
            const int g_start = target_batch_size * gb;
            const int g_end = std::min(target_batch_size * (gb + 1), size_g);
            const int size_gb = g_end - g_start;

            X.resize(size_p, size_gb);
            A.resize(size_p, size_p, size_gb);
            mat_alloc(G_T, spins, size_p, size_gb);

            eval_bfs_per_batch(*xc, X, g_start, g_end);
            SNX_A(*xc, A, Vs, snx_int_thresh, g_start, g_end);

            for(int s = 0; s < spins; s++){
                G_T[s] = contract_A_F(A, *p[s] * X);
                fxc[s] = fxc[s] + X * G_T[s] * (-spin_factor);
            }
        }
        #pragma omp critical
        {
            for(int s = 0; s < spins; s++){
                *(xc->F_XC[s]) = *(xc->F_XC[s]) + fxc[s];
            }
        }
    }
    for(int s = 0; s < spins; s++){
        *xc->F_XC[s] = (*xc->F_XC[s] + transpose(*xc->F_XC[s])) * 0.5; // symmetrize
        xc->E_XC += Tr((*xc->F_XC[s]) * (*p[s]));
    }
	xc->E_XC *= 0.5;

    //auto end_K = std::chrono::high_resolution_clock::now();
    //auto duration_K = std::chrono::duration<double>(end_K - start_K);
    //std::cout << "K time (wall, s) = " << std::setprecision(2) << duration_K.count() << std::setprecision(10) << std::endl;
}

///////////////////////////////////////////////////////////////////
// !!!	 												     !!! //
// !!!    Functions below give PER GRID-POINT F_XC / E_XC    !!! //
// !!!	 												     !!! //
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
// LDA ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

namespace _SLATER{
	inline const double RX = -cbrt(3.0 / M_PI);
	inline const double UX = -cbrt(6.0 / M_PI);
}

void Slater(const XC& xc, const XC_inp& inp, XC_ret& ret){
	using namespace _SLATER;
    if(xc.restricted){
		const double rho   = inp.rho;
		const double rho_3 = cbrt(rho);
		ret.e_XC = RX * (3.0 / 4.0) * rho * rho_3;
		ret.drho_XC = {RX * rho_3};
	}
    else{
		const double rho_a   = inp.rho_a;
		const double rho_b   = inp.rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);	
		ret.e_XC = UX * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		ret.drho_XC = { UX * rho_a_3, UX * rho_b_3 };
	}
}

namespace _VWN5{
	// Paramagnetic (zeta = 0)
	inline const     double A_P  = (1 - log(2)) / (M_PI * M_PI);
	inline constexpr double x0_P = -0.10498;
	inline constexpr double b_P  =  3.72744;
	inline constexpr double c_P  =  12.9352;
	inline constexpr double X0_P =  x0_P * x0_P + b_P * x0_P + c_P;
	inline const     double Q_P  =  sqrt(4 * c_P - b_P * b_P);
	// Ferromagnetic (zeta = 1)
	inline const     double A_F  =  A_P / 2;
	inline constexpr double x0_F = -0.32500;
	inline constexpr double b_F  =  7.06042;
	inline constexpr double c_F  =  18.0578;
	inline constexpr double X0_F =  x0_F * x0_F + b_F * x0_F + c_F;
	inline const     double Q_F  =  sqrt(4 * c_F - b_F * b_F);
	inline const     double ddf0 = 4.0 / (9.0 * (cbrt(2) - 1));
}

void VWN5(const XC& xc, const XC_inp& inp, XC_ret& ret){
	using namespace _SLATER;
	using namespace _VWN5;
	
	if(xc.restricted){
		const double rho = inp.rho;
		const double rho_3 = cbrt(rho);
		const double e_X   = RX * (3.0 / 4.0) * rho * rho_3;
		const double v_X   = RX * rho_3;

		const double x  = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		const double X  = x * x + b_P * x + c_P;

		const double eps_c = A_P * (
			log(x * x / X) + (2 * b_P / Q_P) * (1 - (2 * x0_P + b_P) * x0_P / X0_P) * atan(Q_P / (2 * x + b_P)) - 
			(b_P * x0_P / X0_P) * log((x - x0_P) * (x - x0_P) / X)
		);
		const double v_c = eps_c - A_P * (x / (3 * X)) * (c_P / x - b_P * x0_P / (x - x0_P));

		ret.e_XC = e_X + rho * eps_c;
		ret.drho_XC = { v_X + v_c };
	}
	else{
		const double rho_a   = inp.rho_a;
		const double rho_b   = inp.rho_b;
		const double rho     = rho_a + rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);
		const double e_X     = UX * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		std::vector<double> v = {UX * rho_a_3, UX * rho_b_3};
		
		const double x = sqrt(cbrt(3.0 / (4.0 * M_PI * rho)));
		const double zeta = ( rho_a - rho_b ) / rho;
		const double zeta3 = zeta * zeta * zeta;
		const std::vector<double> dzeta_drho = {2 * rho_b / (rho * rho), -2 * rho_a / (rho * rho)};
		const double f  = f_zeta(zeta);
		const double df = df_zeta(zeta);
		const double alpha = VWN_alpha(x);
		const double dalpha_drho = VWN_dalpha_drho(x, rho);
		const double X_P  = x * x + b_P * x + c_P;
		const double X_F  = x * x + b_F * x + c_F;
		
		const double eps_c_P = A_P * (
			log(x * x / X_P) + (2 * b_P / Q_P) * (1 - (2 * x0_P + b_P) * x0_P / X0_P) * 
			atan(Q_P / (2 * x + b_P)) - (b_P * x0_P / X0_P) * log((x - x0_P) * (x - x0_P) / X_P)
		);
		const double eps_c_F = A_F * (
			log(x * x / X_F) + (2 * b_F / Q_F) * (1 - (2 * x0_F + b_F) * x0_F / X0_F) * 
			atan(Q_F / (2 * x + b_F)) - (b_F * x0_F / X0_F) * log((x - x0_F) * (x - x0_F) / X_F)
		);
		const double e_c = rho * eps_c_P + rho * alpha * (f / ddf0) * (1 - zeta3 * zeta) + rho * (eps_c_F - eps_c_P) * f * zeta3 * zeta;
		
		const double v_c_P = eps_c_P - A_P * (x / (3 * X_P)) * (c_P / x - b_P * x0_P / (x - x0_P));
		const double v_c_F = eps_c_F - A_F * (x / (3 * X_F)) * (c_F / x - b_F * x0_F / (x - x0_F));

		for(int s = 0; s < 2; s++){
			v[s] += v_c_P + (alpha + rho * dalpha_drho) * (f / ddf0) * (1 - zeta3 * zeta) + 
					rho * alpha * ((df/ddf0) * (1 - zeta3 * zeta) - 4 * zeta3 * (f / ddf0)) * dzeta_drho[s] +
					(v_c_F - v_c_P) * f * zeta3 * zeta + 
					rho * (eps_c_F - eps_c_P) * (df * zeta3 * zeta + 4 * zeta3 * f) * dzeta_drho[s];	
		}
		ret.e_XC = e_X + e_c;
		ret.drho_XC = v;
	}
}

namespace _PW92 {
	// Paramagnetic (zeta = 0)
	inline const     double A_P  = (1 - log(2)) / (M_PI * M_PI);
	inline constexpr double a1_P = 0.21370;
	inline constexpr double b1_P = 7.5957;
	inline constexpr double b2_P = 3.5876;
	inline constexpr double b3_P = 1.6382;
	inline constexpr double b4_P = 0.49294;
	// Ferromagnetic (zeta = 1)
	inline const     double A_F  = A_P / 2;
	inline constexpr double a1_F = 0.20548;
	inline constexpr double b1_F = 14.1189;
	inline constexpr double b2_F = 6.1977;
	inline constexpr double b3_F = 3.3662;
	inline constexpr double b4_F = 0.62517;
}

void PW92(const XC& xc, const XC_inp& inp, XC_ret& ret){
	using namespace _SLATER;
	using namespace _PW92;

	if(xc.restricted){
		const double rho = inp.rho;
		const double rho_3 = cbrt(rho);
		const double e_X   = RX * (3.0 / 4.0) * rho * rho_3;
		const double v_X   = RX * rho_3;
		
		const double rs = cbrt(3 / (4 * M_PI * rho));
		const double Q0  = -2 * A_P * (1 + a1_P * rs);
		const double Q1  =  2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs);
		const double Q1p =      A_P * (b1_P / sqrt(rs) + 2 * b2_P + 3 * b3_P * sqrt(rs) + 4 * b4_P * rs);
		
		const double eps_c = -2 * A_P * (1 + a1_P * rs) * log(1 + 
			1 / (2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs)));
		const double v_c = eps_c - (rs / 3) * (-2 * A_P * a1_P * log(1 + 1 / Q1) - Q0 * Q1p / (Q1 * Q1 + Q1));
		
		ret.e_XC = e_X + rho * eps_c;
		ret.drho_XC = { v_X + v_c };
	}
    else{
		const double rho_a   = inp.rho_a;
		const double rho_b   = inp.rho_b;
		const double rho     = rho_a + rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);
		const double e_X     = UX * (3.0 / 4.0) * (rho_a * rho_a_3 + rho_b * rho_b_3);
		std::vector<double> v = {UX * rho_a_3, UX * rho_b_3};
		
		double rs = cbrt(3 / (4 * M_PI * rho));
		const double zeta = (rho_a - rho_b) / rho;
		const double zeta3 = zeta * zeta * zeta;
		const double zeta4 = zeta3 * zeta;
		const double f = f_zeta(zeta);
		const double df = df_zeta(zeta);
		const double alpha = PW92_alpha(rs);
		const double dalpha_drs = PW92_dalpha_drs(rs);
		
		const double eps_P = -2 * A_P * (1 + a1_P * rs) * log(1 + 1 / (2 * A_P * 
			(b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs)));
		const double eps_F = -2 * A_F * (1 + a1_F * rs) * log(1 + 1 / (2 * A_F * 
			(b1_F * sqrt(rs) + b2_F * rs + b3_F * sqrt(intpow(rs, 3)) + b4_F * rs * rs)));
		const double eps = eps_P + alpha * (f / _VWN5::ddf0) * (1 - zeta4) + (eps_F - eps_P) * f * zeta4;
		
		const double Q0_P   = -2 * A_P * (1 + a1_P * rs);
		const double Q1_P   =  2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs);
		const double Q1p_P  =      A_P * (b1_P / sqrt(rs) + 2 * b2_P + 3 * b3_P * sqrt(rs) + 4 * b4_P * rs);
		const double deps_P = -2 * A_P * a1_P * log(1 + 1 / Q1_P) - Q0_P * Q1p_P / (Q1_P * Q1_P + Q1_P);
	
		const double Q0_F   = -2 * A_F * (1 + a1_F * rs);
		const double Q1_F   =  2 * A_F * (b1_F * sqrt(rs) + b2_F * rs + b3_F * sqrt(intpow(rs, 3)) + b4_F * rs * rs);
		const double Q1p_F  =      A_F * (b1_F / sqrt(rs) + 2 * b2_F + 3 * b3_F * sqrt(rs) + 4 * b4_F * rs);
		const double deps_F = -2 * A_F * a1_F * log(1 + 1 / Q1_F) - Q0_F * Q1p_F / (Q1_F * Q1_F + Q1_F);
		
		const double deps_dr = deps_P * (1 - f * zeta4) + deps_F * f * zeta4 + dalpha_drs * (f / _VWN5::ddf0) * (1 - zeta4);
		const double deps_dz = 4 * zeta3 * f * (eps_F - eps_P - alpha / _VWN5::ddf0) + 
			df * (zeta4 * (eps_F - eps_P) + (1 - zeta4) * alpha / _VWN5::ddf0);
		
		v[0] += eps - (rs / 3) * deps_dr - (zeta - 1) * deps_dz;
		v[1] += eps - (rs / 3) * deps_dr - (zeta + 1) * deps_dz;

		ret.e_XC = e_X + rho * eps;
		ret.drho_XC = v;
	}
}

///////////////////////////////////////////////////////////////////
// GGA ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

namespace _PBE{
	inline constexpr double beta = 0.066725;
	inline constexpr double kappa = 0.804;
	inline constexpr double mu = beta * (M_PI * M_PI / 3);
	inline const     double gamma = (1 - log(2)) / (M_PI * M_PI);
}

void PBE_X(const XC& xc, const XC_inp& inp, XC_ret& ret){
	using namespace _SLATER;
	using namespace _PBE;
	
	if(xc.restricted){
		// Slater Exchange
		const double rho = inp.rho;
		const double rho_3 = cbrt(rho);
		const double e_LDA = RX * (3.0 / 4.0) * rho * rho_3;
		const double v_LDA = RX * rho_3;

		// Enhancement Factor
		const double grho2 = (
			inp.gradient_rho[0] * inp.gradient_rho[0] + 
			inp.gradient_rho[1] * inp.gradient_rho[1] + 
			inp.gradient_rho[2] * inp.gradient_rho[2] 
		);
		const double kF = cbrt(3 * M_PI * M_PI * rho);
		const double s2 = grho2 / (4 * kF * kF * rho * rho);
		const double ds2_drho = -8.0 * s2 / (3.0 * rho);
		const double ds2_dgrho2 = s2 / grho2;	
		const double FX_d = 1 + mu * s2 / kappa;	
		const double FX = 1 + kappa - kappa / FX_d;
		const double dFX_ds2 = mu / (FX_d * FX_d);
		const double dFX_drho = dFX_ds2 * ds2_drho; 
		const double dFX_dgrho2 = dFX_ds2 * ds2_dgrho2;
 
		const double de_drho = v_LDA * FX + e_LDA * dFX_drho;
		const double de_dgrho2 = e_LDA * dFX_dgrho2;

		ret.e_XC = e_LDA * FX;
		ret.drho_XC = {de_drho};
		ret.dgamma_XC = {de_dgrho2};
	}
    else{
		// Slater Exchange
		const double rho_a = 2 * inp.rho_a;
		const double rho_b = 2 * inp.rho_b;
		const double rho_a_3 = cbrt(rho_a);
		const double rho_b_3 = cbrt(rho_b);
		const double e_LDA_a = (3.0 / 4.0) * RX * rho_a * rho_a_3;
		const double e_LDA_b = (3.0 / 4.0) * RX * rho_b * rho_b_3;
		const double v_LDA_a = RX * rho_a_3;	// v_LDA_X spin scaling already accounted for!
		const double v_LDA_b = RX * rho_b_3;	// v_LDA_X spin scaling already accounted for!
		
		// Enhancement Factor
		const double grho2_a = 4 * (
			inp.gradient_rho_a[0] * inp.gradient_rho_a[0] + 
			inp.gradient_rho_a[1] * inp.gradient_rho_a[1] + 
			inp.gradient_rho_a[2] * inp.gradient_rho_a[2]
		); 
		const double grho2_b = 4 * (
			inp.gradient_rho_b[0] * inp.gradient_rho_b[0] + 
			inp.gradient_rho_b[1] * inp.gradient_rho_b[1] + 
			inp.gradient_rho_b[2] * inp.gradient_rho_b[2]
		); 
		const double kF_a = cbrt(3 * M_PI * M_PI * rho_a);
		const double kF_b = cbrt(3 * M_PI * M_PI * rho_b);
		const double s2_a = grho2_a / (4 * kF_a * kF_a * rho_a * rho_a);
		const double s2_b = grho2_b / (4 * kF_b * kF_b * rho_b * rho_b);
		const double ds2_drho_a = -8.0 * s2_a / (3.0 * rho_a);
		const double ds2_drho_b = -8.0 * s2_b / (3.0 * rho_b);
		const double ds2_dgrho2_a = s2_a / grho2_a;	
		const double ds2_dgrho2_b = s2_b / grho2_b;	
		const double FX_d_a = 1 + mu * s2_a / kappa;	
		const double FX_d_b = 1 + mu * s2_b / kappa;	
		const double FX_a = 1 + kappa - kappa / FX_d_a;
		const double FX_b = 1 + kappa - kappa / FX_d_b;
		const double dFX_ds2_a = mu / (FX_d_a * FX_d_a);
		const double dFX_ds2_b = mu / (FX_d_b * FX_d_b);
		const double dFX_drho_a = dFX_ds2_a * ds2_drho_a;
		const double dFX_drho_b = dFX_ds2_b * ds2_drho_b;
		const double dFX_dgrho2_a = dFX_ds2_a * ds2_dgrho2_a;
		const double dFX_dgrho2_b = dFX_ds2_b * ds2_dgrho2_b;
		
		const double de_drho_a = v_LDA_a * FX_a + e_LDA_a * dFX_drho_a;	// All divided by
		const double de_drho_b = v_LDA_b * FX_b + e_LDA_b * dFX_drho_b;	// 2 to account
		const double de_dgrho2_a = 2 * e_LDA_a * dFX_dgrho2_a;			// for spin scaling
		const double de_dgrho2_b = 2 * e_LDA_b * dFX_dgrho2_b;			// of e_X and v_X

		ret.e_XC = 0.5 * (e_LDA_a * FX_a + e_LDA_b * FX_b);
		ret.drho_XC = {de_drho_a, de_drho_b};
		ret.dgamma_XC = {de_dgrho2_a, de_dgrho2_b};
	}
}

void PBE(const XC& xc, const XC_inp& inp, XC_ret& ret){
	using namespace _SLATER;
	using namespace _PW92;
	using namespace _PBE;
	assert(xc.restricted); // only restricted PBE for now!
	
    // Slater Exchange
    const double rho = inp.rho;
    const double rho_3 = cbrt(rho);
    const double e_LDA_X = RX * (3.0 / 4.0) * rho * rho_3;
    const double v_LDA_X = RX * rho_3;

    // Enhancement Factor
    const double grho2 = (
        inp.gradient_rho[0] * inp.gradient_rho[0] + 
        inp.gradient_rho[1] * inp.gradient_rho[1] + 
        inp.gradient_rho[2] * inp.gradient_rho[2] 
    );
    const double kF = cbrt(3 * M_PI * M_PI * rho);
    const double s2 = grho2 / (4 * kF * kF * rho * rho);
    const double ds2_drho = -8.0 * s2 / (3.0 * rho);
    const double ds2_dgrho2 = s2 / grho2;	
    const double FX_d = 1 + mu * s2 / kappa;	
    const double FX = 1 + kappa - kappa / FX_d;
    const double dFX_ds2 = mu / (FX_d * FX_d);
    const double dFX_drho = dFX_ds2 * ds2_drho; 
    const double dFX_dgrho2 = dFX_ds2 * ds2_dgrho2;

    const double de_X_drho = v_LDA_X * FX + e_LDA_X * dFX_drho;
    const double de_X_dgrho2 = e_LDA_X * dFX_dgrho2;

    ret.e_XC = e_LDA_X * FX;
    ret.drho_XC = {de_X_drho};
    ret.dgamma_XC = {de_X_dgrho2};

    // PW92 correlation
    const double rs = cbrt(3 / (4 * M_PI * rho));
    const double Q0  = -2 * A_P * (1 + a1_P * rs);
    const double Q1  =  2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs);
    const double Q1p =      A_P * (b1_P / sqrt(rs) + 2 * b2_P + 3 * b3_P * sqrt(rs) + 4 * b4_P * rs);
    
    const double eps_c_LDA = -2 * A_P * (1 + a1_P * rs) * log(1 + 
        1 / (2 * A_P * (b1_P * sqrt(rs) + b2_P * rs + b3_P * sqrt(intpow(rs, 3)) + b4_P * rs * rs)));
    const double deps_c_LDA_dn = -(rs / 3) * (-2 * A_P * a1_P * log(1 + 1 / Q1) - Q0 * Q1p / (Q1 * Q1 + Q1)) / rho;
    const double v_c_LDA = eps_c_LDA + rho * deps_c_LDA_dn;

    // PBE correlation correction
    const double ks = sqrt(4 * kF / M_PI);
    const double t2 = grho2 / (4 * ks * ks * rho * rho);
    double A_PBE = (beta / gamma) / (exp(-eps_c_LDA / gamma) - 1);
    A_PBE = (A_PBE > 1e10 ? 1e10 : A_PBE);
    const double dnm = 1 + A_PBE * t2 + A_PBE * A_PBE * t2 * t2;
    const double Q = 1 + (beta / gamma) * t2 * (1 + A_PBE * t2) / dnm;
    const double H = gamma * log(Q);

    const double dH_dt2 = (beta / Q) * (1 + 2 * A_PBE * t2) / (dnm * dnm);
    const double dt2_dn = -(7.0 / 3.0) * t2 / rho;
    const double dH_dA_PBE = -(beta / Q) * A_PBE * t2 * t2 * t2 * (2 + A_PBE * t2) / (dnm * dnm);
    const double dA_PBE_deps = A_PBE * (A_PBE + beta / gamma) / beta;
    const double dt2_dgrho2 = t2 / grho2;

    ret.e_XC += rho * (eps_c_LDA + H);
    ret.drho_XC[0] += v_c_LDA + H + rho * (dH_dt2 * dt2_dn + dH_dA_PBE * dA_PBE_deps * deps_c_LDA_dn);
    ret.dgamma_XC[0] += rho * dH_dt2 * dt2_dgrho2;
}

///////////////////////////////////////////////////////////////////
// MGGA ///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

namespace _B97M_V{
	inline const double c_tau_ueg = 3.0 * cbrt(36 * M_PI) * M_PI / 5.0;

	inline constexpr double c_x00 = 1.000;
	inline constexpr double c_x10 = 0.416;
	inline constexpr double c_x01 = 1.308;
	inline constexpr double c_x11 = 3.070;
	inline constexpr double c_x02 = 1.901;

	inline constexpr double c_css00 =  1.000;
	inline constexpr double c_css10 = -5.668;
	inline constexpr double c_css02 = -1.855;
	inline constexpr double c_css32 = -20.497;
	inline constexpr double c_css42 = -20.364;

	inline constexpr double c_cos00 =  1.000;
	inline constexpr double c_cos10 =  2.535;
	inline constexpr double c_cos01 =  1.573;
	inline constexpr double c_cos32 = -6.427;
	inline constexpr double c_cos03 = -6.298;

	inline constexpr double gamma_x   = 0.004;
	inline constexpr double gamma_css = 0.2;
	inline constexpr double gamma_cos = 0.006;

    inline constexpr double b_VV10 = 6.00;
    inline constexpr double C_VV10 = 0.01;
}

void B97M_V(const XC& xc, const XC_inp& inp, XC_ret& ret){
	using namespace _SLATER;
	using namespace _B97M_V;
	assert(!xc.restricted);
	// Separation of same/opp spin correlation 
	// requires unrestricted calculation!
    // E = e_X + e_css + e_cos + e_VV10
    constexpr double DIV_0_GUARD = 1e-20;

    const double rho   = inp.rho;
    const double rho_a = inp.rho_a;
    const double rho_b = inp.rho_b;
    const double rho_a_div = (rho_a > DIV_0_GUARD ? rho_a : DIV_0_GUARD);
    const double rho_b_div = (rho_b > DIV_0_GUARD ? rho_b : DIV_0_GUARD);

    //const std::vector<double>& grho = inp.gradient_rho;
    const std::vector<double>& grho_a = inp.gradient_rho_a;
    const std::vector<double>& grho_b = inp.gradient_rho_b;
    //const double grho2 = grho[0] * grho[0] + grho[1] * grho[1] + grho[2] * grho[2];
    const double grho2_a = grho_a[0] * grho_a[0] + grho_a[1] * grho_a[1] + grho_a[2] * grho_a[2];
    const double grho2_b = grho_b[0] * grho_b[0] + grho_b[1] * grho_b[1] + grho_b[2] * grho_b[2];
    const double grho2_a_div = (grho2_a > DIV_0_GUARD ? grho2_a : DIV_0_GUARD);
    const double grho2_b_div = (grho2_b > DIV_0_GUARD ? grho2_b : DIV_0_GUARD);

    const double tau_a = inp.ke_density_a;
    const double tau_b = inp.ke_density_b;
    const double tau_a_div = (tau_a > DIV_0_GUARD ? tau_a : DIV_0_GUARD);
    const double tau_b_div = (tau_b > DIV_0_GUARD ? tau_b : DIV_0_GUARD);
    const double tau_ueg_a = c_tau_ueg * rho_a * cbrt(rho_a * rho_a);
    const double tau_ueg_b = c_tau_ueg * rho_b * cbrt(rho_b * rho_b);
    const double t_a = tau_ueg_a / tau_a_div;
    const double t_b = tau_ueg_b / tau_b_div;

    const double dt_a_drho_a = (5.0 / 3.0) * t_a / rho_a_div;
    const double dt_b_drho_b = (5.0 / 3.0) * t_b / rho_b_div;
    const double dt_a_dtau_a = -t_a / tau_a_div;
    const double dt_b_dtau_b = -t_b / tau_b_div;

    // e_X
    const double e_X_ueg_a = UX * (3.0 / 4.0) * (rho_a * cbrt(rho_a));
    const double e_X_ueg_b = UX * (3.0 / 4.0) * (rho_b * cbrt(rho_b));
    const double s2_a = grho2_a / (rho_a_div * rho_a_div * cbrt(rho_a_div * rho_a_div));
    const double s2_b = grho2_b / (rho_b_div * rho_b_div * cbrt(rho_b_div * rho_b_div));
    const double wx_a = (t_a - 1) / (t_a + 1);
    const double wx_b = (t_b - 1) / (t_b + 1);
    const double ux_a = gamma_x * s2_a / (1 + gamma_x * s2_a);
    const double ux_b = gamma_x * s2_b / (1 + gamma_x * s2_b);

    const double gx_a = c_x00 + (c_x10 * wx_a) + (c_x01 * ux_a) + (c_x11 * wx_a * ux_a) + (c_x02 * ux_a * ux_a);
    const double gx_b = c_x00 + (c_x10 * wx_b) + (c_x01 * ux_b) + (c_x11 * wx_b * ux_b) + (c_x02 * ux_b * ux_b);

    const double e_X = e_X_ueg_a * gx_a + e_X_ueg_b * gx_b;

    // e_X derivatives
    const double de_X_ueg_a_drho_a = UX * cbrt(rho_a);
    const double de_X_ueg_b_drho_b = UX * cbrt(rho_b);
    const double ds2_a_drho_a = -(8.0 / 3.0) * s2_a / rho_a_div;
    const double ds2_b_drho_b = -(8.0 / 3.0) * s2_b / rho_b_div;
    const double ds2_a_dgrho2_a = s2_a / grho2_a_div;
    const double ds2_b_dgrho2_b = s2_b / grho2_b_div;
    const double dwx_a_dt_a = 2.0 / ((t_a + 1) * (t_a + 1));
    const double dwx_b_dt_b = 2.0 / ((t_b + 1) * (t_b + 1));
    const double dux_a_ds2_a = gamma_x / ((1 + gamma_x * s2_a) * (1 + gamma_x * s2_a));
    const double dux_b_ds2_b = gamma_x / ((1 + gamma_x * s2_b) * (1 + gamma_x * s2_b));
    const double dgx_a_dwx_a = c_x10 + (c_x11 * ux_a);	
    const double dgx_b_dwx_b = c_x10 + (c_x11 * ux_b);
    const double dgx_a_dux_a = c_x01 + (c_x11 * wx_a) + (2 * c_x02 * ux_a);
    const double dgx_b_dux_b = c_x01 + (c_x11 * wx_b) + (2 * c_x02 * ux_b);

    const double de_X_drho_a = de_X_ueg_a_drho_a * gx_a + 
        e_X_ueg_a * (dgx_a_dwx_a * dwx_a_dt_a * dt_a_drho_a + dgx_a_dux_a * dux_a_ds2_a * ds2_a_drho_a);
    const double de_X_drho_b = de_X_ueg_b_drho_b * gx_b + 
        e_X_ueg_b * (dgx_b_dwx_b * dwx_b_dt_b * dt_b_drho_b + dgx_b_dux_b * dux_b_ds2_b * ds2_b_drho_b);

    const double de_X_dgrho2_a = e_X_ueg_a * dgx_a_dux_a * dux_a_ds2_a * ds2_a_dgrho2_a;
    const double de_X_dgrho2_b = e_X_ueg_b * dgx_b_dux_b * dux_b_ds2_b * ds2_b_dgrho2_b;

    const double de_X_dtau_a = e_X_ueg_a * dgx_a_dwx_a * dwx_a_dt_a * dt_a_dtau_a;
    const double de_X_dtau_b = e_X_ueg_b * dgx_b_dwx_b * dwx_b_dt_b * dt_b_dtau_b;

    // e_css
    const double eps_pw92_aa = eps_c_pw92(rho_a, 0.0);
    const double eps_pw92_bb = eps_c_pw92(0.0, rho_b);
    const double wc_aa = wx_a;
    const double wc_bb = wx_b;
    const double uc_aa = gamma_css * s2_a / (1 + gamma_css * s2_a);
    const double uc_bb = gamma_css * s2_b / (1 + gamma_css * s2_b);

    const double gcss_aa = c_css00 + (c_css10 * wc_aa) + (c_css02 * uc_aa * uc_aa) + (c_css32 * intpow(wc_aa, 3) * uc_aa * uc_aa) + 
        (c_css42 * intpow(wc_aa, 4) * uc_aa * uc_aa);
    const double gcss_bb = c_css00 + (c_css10 * wc_bb) + (c_css02 * uc_bb * uc_bb) + (c_css32 * intpow(wc_bb, 3) * uc_bb * uc_bb) + 
        (c_css42 * intpow(wc_bb, 4) * uc_bb * uc_bb);

    const double e_css = (rho_a * eps_pw92_aa * gcss_aa) + (rho_b * eps_pw92_bb * gcss_bb);

    // e_css_derivatives
    const double deps_pw92_aa_drho_a = deps_c_dns_pw92(rho_a, 0.0, 0);
    const double deps_pw92_bb_drho_b = deps_c_dns_pw92(0.0, rho_b, 1);
    const double dwc_aa_dt_a = dwx_a_dt_a;
    const double dwc_bb_dt_b = dwx_b_dt_b;
    const double duc_aa_ds2_a = gamma_css / ((1 + gamma_css * s2_a) * (1 + gamma_css * s2_a));
    const double duc_bb_ds2_b = gamma_css / ((1 + gamma_css * s2_b) * (1 + gamma_css * s2_b));
    
    const double dgcss_aa_dwc_aa = c_css10 + (3 * c_css32 * wc_aa * wc_aa * uc_aa * uc_aa) + 
        (4 * c_css42 * intpow(wc_aa, 3) * uc_aa * uc_aa);
    const double dgcss_bb_dwc_bb = c_css10 + (3 * c_css32 * wc_bb * wc_bb * uc_bb * uc_bb) + 
        (4 * c_css42 * intpow(wc_bb, 3) * uc_bb * uc_bb);

    const double dgcss_aa_duc_aa = (2 * c_css02 * uc_aa) + (2 * c_css32 * intpow(wc_aa, 3) * uc_aa) + 
        (2 * c_css42 * intpow(wc_aa, 4) * uc_aa); 
    const double dgcss_bb_duc_bb = (2 * c_css02 * uc_bb) + (2 * c_css32 * intpow(wc_bb, 3) * uc_bb) + 
        (2 * c_css42 * intpow(wc_bb, 4) * uc_bb);

    const double de_css_drho_a = (eps_pw92_aa + rho_a * deps_pw92_aa_drho_a) * gcss_aa + 
        rho_a * eps_pw92_aa * (dgcss_aa_dwc_aa * dwc_aa_dt_a * dt_a_drho_a + dgcss_aa_duc_aa * duc_aa_ds2_a * ds2_a_drho_a);
    const double de_css_drho_b = (eps_pw92_bb + rho_b * deps_pw92_bb_drho_b) * gcss_bb + 
        rho_b * eps_pw92_bb * (dgcss_bb_dwc_bb * dwc_bb_dt_b * dt_b_drho_b + dgcss_bb_duc_bb * duc_bb_ds2_b * ds2_b_drho_b);

    const double de_css_dgrho2_a = rho_a * eps_pw92_aa * (dgcss_aa_duc_aa * duc_aa_ds2_a * ds2_a_dgrho2_a);
    const double de_css_dgrho2_b = rho_b * eps_pw92_bb * (dgcss_bb_duc_bb * duc_bb_ds2_b * ds2_b_dgrho2_b);

    const double de_css_dtau_a = rho_a * eps_pw92_aa * (dgcss_aa_dwc_aa * dwc_aa_dt_a * dt_a_dtau_a);
    const double de_css_dtau_b = rho_b * eps_pw92_bb * (dgcss_bb_dwc_bb * dwc_bb_dt_b * dt_b_dtau_b);

    // e_cos
    const double t_ab = 0.5 * (t_a + t_b);
    const double s2ab = 0.5 * (s2_a + s2_b);
    const double e_pw92_ab = rho * eps_c_pw92(rho_a, rho_b) - rho_a * eps_pw92_aa - rho_b * eps_pw92_bb;
    const double wc_ab = (t_ab - 1) / (t_ab + 1);
    const double uc_ab = gamma_cos * s2ab / (1 + gamma_cos * s2ab);

    const double gcos = c_cos00 + (c_cos10 * wc_ab) + (c_cos01 * uc_ab) + (c_cos32 * intpow(wc_ab, 3) * uc_ab * uc_ab) + 
        (c_cos03 * intpow(uc_ab, 3));

    const double e_cos = e_pw92_ab * gcos;

    // e_cos derivatives
    const double de_pw92_ab_drho_a = (eps_c_pw92(rho_a, rho_b) + rho * deps_c_dns_pw92(rho_a, rho_b, 0)) - 
        (eps_pw92_aa + rho_a * deps_pw92_aa_drho_a);
    const double de_pw92_ab_drho_b = (eps_c_pw92(rho_a, rho_b) + rho * deps_c_dns_pw92(rho_a, rho_b, 1)) - 
        (eps_pw92_bb + rho_b * deps_pw92_bb_drho_b);
    const double dwc_ab_dt_ab = 2.0 / ((t_ab + 1) * (t_ab + 1));
    const double duc_ab_ds2_ab = gamma_cos / ((1 + gamma_cos * s2ab) * (1 + gamma_cos * s2ab));
    
    const double dgcos_dwc_ab = c_cos10 + (3 * c_cos32 * wc_ab * wc_ab * uc_ab * uc_ab);
    const double dgcos_duc_ab = c_cos01 + (2 * c_cos32 * wc_ab * wc_ab * wc_ab * uc_ab) + (3 * c_cos03 * uc_ab * uc_ab);

    const double de_cos_drho_a = de_pw92_ab_drho_a * gcos + e_pw92_ab * 
        (dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_a_drho_a + dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_a_drho_a);
    const double de_cos_drho_b = de_pw92_ab_drho_b * gcos + e_pw92_ab * 
        (dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_b_drho_b + dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_b_drho_b);

    const double de_cos_dgrho2_a = e_pw92_ab * (dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_a_dgrho2_a);
    const double de_cos_dgrho2_b = e_pw92_ab * (dgcos_duc_ab * duc_ab_ds2_ab * 0.5 * ds2_b_dgrho2_b);

    const double de_cos_dtau_a = e_pw92_ab * (dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_a_dtau_a);
    const double de_cos_dtau_b = e_pw92_ab * (dgcos_dwc_ab * dwc_ab_dt_ab * 0.5 * dt_b_dtau_b);

    ret.e_XC = e_X + e_css + e_cos;
    ret.drho_XC = {
        (de_X_drho_a + de_css_drho_a + de_cos_drho_a),
        (de_X_drho_b + de_css_drho_b + de_cos_drho_b)
    };
    ret.dgamma_XC = {
        (de_X_dgrho2_a + de_css_dgrho2_a + de_cos_dgrho2_a), 
        (de_X_dgrho2_b + de_css_dgrho2_b + de_cos_dgrho2_b)
    };
    ret.dtau_XC = {
        (de_X_dtau_a + de_css_dtau_a + de_cos_dtau_a), 
        (de_X_dtau_b + de_css_dtau_b + de_cos_dtau_b)
    };
    // VV10 added later, see VV10 below
}

void eval_VV10_properties(const XC& xc, Matrix& phi_buf, Matrix& gpx_buf, Matrix& gpy_buf, Matrix& gpz_buf, 
    std::vector<double>& rhos, std::vector<std::vector<double>>& grds)
{
    const int size_p = xc.mol->AOs.size();
    const int size_g = xc.g->num_gridpoints;
    const int size_r = rhos.size();
    const int size_d = grds.size();

    assert((phi_buf.rows==size_g) && (phi_buf.cols==size_p));
    assert((gpx_buf.rows==size_g) && (gpx_buf.cols==size_p));
    assert((gpy_buf.rows==size_g) && (gpy_buf.cols==size_p));
    assert((gpz_buf.rows==size_g) && (gpz_buf.cols==size_p));
    assert((size_r==size_g) && (size_d==size_g));

    const std::vector<GF> bfs = xc.mol->AOs;
    const std::vector<double>& gx = xc.g->x;
    const std::vector<double>& gy = xc.g->y;
    const std::vector<double>& gz = xc.g->z;

    const int spinidx = (xc.restricted ? 0 : 2);
    const Matrix& p = *(xc.P[spinidx]);
    
    #pragma omp parallel
    {
        std::vector<double> tmp_grd(3);
        #pragma omp for
        for(int g = 0; g < size_g; g++){
            for(int j = 0; j < size_p; j++){
                phi_buf(g, j) = bfs[j].evaluate(gx[g], gy[g], gz[g]);

                tmp_grd = bfs[j].evaluate_gradient(gx[g], gy[g], gz[g]);
                gpx_buf(g, j) = tmp_grd[0];
                gpy_buf(g, j) = tmp_grd[1];
                gpz_buf(g, j) = tmp_grd[2];
            }
            rhos[g] = density(phi_buf.getRow(g), p);
            grds[g] = density_gradient(phi_buf.getRow(g), gpx_buf.getRow(g), gpy_buf.getRow(g), gpz_buf.getRow(g), p);
        }
    }
}

void VV10(XC* xc){
	constexpr double DIV_0_GUARD = 1e-20;
    const double b = xc->nlc_params[0];
    const double C = xc->nlc_params[1];
    const double beta = sqrt(sqrt(27.0) / b) / 32.0 / b;
    
    const int size_p = xc->mol->AOs.size();
    const int size_g = xc->g->num_gridpoints;

    Matrix phi_buf(size_g, size_p);
    Matrix gpx_buf(size_g, size_p);
    Matrix gpy_buf(size_g, size_p);
    Matrix gpz_buf(size_g, size_p);
    std::vector<double> rho(size_g);
    std::vector<std::vector<double>> grd(size_g);

    eval_VV10_properties(*xc, phi_buf, gpx_buf, gpy_buf, gpz_buf, rho, grd);

    const std::vector<double>& gx = xc->g->x;
    const std::vector<double>& gy = xc->g->y;
    const std::vector<double>& gz = xc->g->z;
    const std::vector<double>& gw = xc->g->w;
   
    const int spins = (xc->restricted ? 1 : 2); 
    double E_XC  = 0.0;
    #pragma omp parallel
    {
        std::vector<double> ref_grd(3);
        std::vector<double> grd_gpt(3);
        Matrix F_XC(size_p, size_p);
        #pragma omp for reduction(+:E_XC)
        for(int g = 0; g < size_g * size_g; g++){
            const int rg = g / size_g;
            const int sg = g % size_g;

            const double ref_rho = rho[rg];
            const double ref_rho_div = (ref_rho > DIV_0_GUARD ? ref_rho : DIV_0_GUARD);
            const double rho_gpt = rho[sg];
            const double rho_gpt_div = (rho_gpt > DIV_0_GUARD ? rho_gpt : DIV_0_GUARD);
            if((ref_rho < 1e-20) || (rho_gpt < 1e-20)){continue;}

            ref_grd = grd[rg];
            const double ref_grho2 = (
                ref_grd[0] * ref_grd[0] +
                ref_grd[1] * ref_grd[1] +
                ref_grd[2] * ref_grd[2]
            );
            const double ref_grho2_div = (ref_grho2 > DIV_0_GUARD ? ref_grho2 : DIV_0_GUARD);
            grd_gpt = grd[sg];
            const double grho2_gpt = (
                grd_gpt[0] * grd_gpt[0] +
                grd_gpt[1] * grd_gpt[1] +
                grd_gpt[2] * grd_gpt[2]
            );	
            const double ref_omega_p2 = 4.0 * M_PI * ref_rho;
            const double ref_omega_g2 = C * ref_grho2 * ref_grho2 / (ref_rho_div * ref_rho_div * ref_rho_div * ref_rho_div);
            const double ref_omega_0 = sqrt(ref_omega_g2 + ref_omega_p2 / 3.0);
            const double ref_kappa = b * (3.0 * M_PI / 2.0) * sqrt(cbrt(ref_rho / (9.0 * M_PI)));
            
            const double omega_p2 = 4.0 * M_PI * rho_gpt;
            const double omega_g2 = C * grho2_gpt * grho2_gpt / (rho_gpt_div * rho_gpt_div * rho_gpt_div * rho_gpt_div);
            const double omega_0 = sqrt(omega_g2 + omega_p2 / 3.0);
            const double kappa = b * (3.0 * M_PI / 2.0) * sqrt(cbrt(rho_gpt / (9.0 * M_PI)));

            const double R2 = intpow(gx[rg] - gx[sg], 2) + intpow(gy[rg] - gy[sg], 2) + intpow(gz[rg] - gz[sg], 2);
            const double ref_g   = ref_omega_0 * R2 + ref_kappa;
            const double g_prime = omega_0 * R2 + kappa;

            const double PHI = -3.0 / (2.0 * ref_g * g_prime * (ref_g + g_prime));
            E_XC += gw[rg] * gw[sg] * 0.5 * ref_rho * rho_gpt * PHI;
            
            const double ref_dkappa_drho = ref_kappa / (6.0 * ref_rho_div);
            const double ref_domega_0_drho = (2.0 / ref_omega_0) * (M_PI / 3.0 - ref_omega_g2 / ref_rho_div);
            const double g_gprime_term = 1.0 / ref_g + 1.0 / (ref_g + g_prime);
            const double FXC_LDA_term = gw[rg] * gw[sg] * rho_gpt * PHI * (
                    1.0 - ref_rho * g_gprime_term * (ref_dkappa_drho + R2 * ref_domega_0_drho)
                ); 
            const double ref_domega_0_dgamma = ref_omega_g2 / (ref_grho2_div * ref_omega_0);
            const double FXC_GGA_term = -gw[rg] * gw[sg] * ref_rho * rho_gpt * ref_domega_0_dgamma * R2 * PHI * g_gprime_term;

            for(int mu = 0; mu < size_p; mu++){
                F_XC(mu, mu) += phi_buf(rg, mu) * FXC_LDA_term * phi_buf(rg, mu) + 
                    2 * FXC_GGA_term * (
		                phi_buf(rg, mu) * (ref_grd[0] * gpx_buf(rg, mu) + ref_grd[1] * gpy_buf(rg, mu) + ref_grd[2] * gpz_buf(rg, mu)) + 
		                phi_buf(rg, mu) * (ref_grd[0] * gpx_buf(rg, mu) + ref_grd[1] * gpy_buf(rg, mu) + ref_grd[2] * gpz_buf(rg, mu))
                    );
                for(int nu = 0; nu < mu; nu++){
                    F_XC(mu, nu) += phi_buf(rg, mu) * FXC_LDA_term * phi_buf(rg, nu) + 
                        2 * FXC_GGA_term * (
		                    phi_buf(rg, nu) * (ref_grd[0] * gpx_buf(rg, mu) + ref_grd[1] * gpy_buf(rg, mu) + ref_grd[2] * gpz_buf(rg, mu)) + 
		                    phi_buf(rg, mu) * (ref_grd[0] * gpx_buf(rg, nu) + ref_grd[1] * gpy_buf(rg, nu) + ref_grd[2] * gpz_buf(rg, nu))
                        );
                    F_XC(nu, mu) = F_XC(mu, nu);
                }
            }
        }
        #pragma omp critical
        {
            for(int s = 0; s < spins; s++){
                *(xc->F_XC[s]) = *(xc->F_XC[s]) + F_XC;
            }
        }
    }
    for(int s = 0; s < spins; s++){
        *(xc->F_XC[s]) = *(xc->F_XC[s]) + *(xc->overlap) * beta;
    }
    E_XC += beta * xc->mol->Nelec;
    xc->E_XC += E_XC;
}

