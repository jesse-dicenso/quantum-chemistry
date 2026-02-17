#include "snx_helper.hpp"
#include "func.hpp"

// Build 3 index tensor A_{\nu \lambda g} for batch of gridpoints
void SNX_A(const XC& xc, Tensor3& A, const Matrix& Vs, double int_thresh, int g_start, int g_end){
    const int size_bf = xc.mol->AOs.size();
    assert((A.dim1==size_bf) && (A.dim2==size_bf) && (A.dim3==(g_end-g_start)));
	const std::vector<GF>& bfs = xc.mol->AOs;
    const std::vector<double>& x = xc.g->x;
    const std::vector<double>& y = xc.g->y;
    const std::vector<double>& z = xc.g->z;
    const std::vector<double>& w = xc.g->w;
    std::vector<double> xyz_g(3);
	for(int nu = 0; nu < size_bf; nu++){
		for(int ld = 0; ld < size_bf; ld++){
            if(Vs(nu, ld) < int_thresh){continue;}
            for(int g = g_start; g < g_end; g++){
                xyz_g[0] = x[g];
                xyz_g[1] = y[g];
                xyz_g[2] = z[g];
		    	A(nu, ld, g-g_start) = w[g] * V(bfs[nu], bfs[ld], xyz_g); // w[g] included in A
		    }
	    }
    }
}

// Contract A and F to give G^T, assumes g is fast idx of A
Matrix contract_A_F(const Tensor3& A, const Matrix& F){
    const int size_bf = F.rows;
    const int size_bg = A.dim3;
    assert((A.dim1 == size_bf) && (A.dim2 == size_bf));
	Matrix G_T(size_bg, size_bf); // G is size_bf * size_bg and this will be G^T
	for(int nu = 0; nu < size_bf; nu++){
		for(int ld = 0; ld < size_bf; ld++){
	        for(int g = 0; g < size_bg; g++){
				G_T(g, nu) += A(nu, ld, g) * F(ld, g);
			}
		}
	}
	return G_T;
}

Matrix V_screen(const std::vector<GF>& bfs){
    const int size_bf = bfs.size();
    Matrix result(size_bf, size_bf);
    for(int mu = 0; mu < size_bf; mu++){
        result(mu, mu) = V_screen_cGTOs(bfs[mu], bfs[mu]);
        //std::cout << "V" << mu << "," << mu << " = " << result(mu,mu) << std::endl;
        for(int nu = 0; nu < mu; nu++){
            result(mu, nu) = V_screen_cGTOs(bfs[mu], bfs[nu]);
            //std::cout << "V" << mu << "," << nu << " = " << result(mu,nu) << std::endl;
            result(nu, mu) = result(mu, nu);
        }
    }
    return result;
}

double V_screen_cGTOs(const GF& g1, const GF& g2){
    const int size_g1 = g1.exps.size();
    const int size_g2 = g2.exps.size();
    
    const int l1  = g1.shell[0] + g1.shell[1] + g1.shell[2];
    const int l1p = (l1%2==0 ? l1 : l1+1);
    const int l1m = (l1%2==0 ? l1 : l1-1);

    const int l2 = g2.shell[0] + g2.shell[1] + g2.shell[2];
    const int l2p = (l2%2==0 ? l2 : l2+1);
    const int l2m = (l2%2==0 ? l2 : l2-1);

    const double AB = sqrt(intpow(g1.xyz[0]-g2.xyz[0], 2) + intpow(g1.xyz[1]-g2.xyz[1], 2) + intpow(g1.xyz[2]-g2.xyz[2], 2));

    double sum = 0.0;
    for(int p = 0; p < size_g1; p++){
        for(int q = 0; q < size_g2; q++){
            sum += g1.N[p] * g2.N[q] * fabs(g1.d[p] * g2.d[q]) * sqrt(
                    V_screen_primitives(l1p, l2m, g1.exps[p], g2.exps[q], AB) * 
                    V_screen_primitives(l1m, l2p, g1.exps[p], g2.exps[q], AB)
                );
        }
    }
    return sum;
}

// Eqn A(22) in 10.1063/1.5048491
double V_screen_primitives(int r, int s, double a, double b, double QAB){
    auto B_Ylm = [](double exp, int l){
        const double numr = sqrt(intpow(2 * exp, 2 * l + 3));
        const double denm = tgamma(l / 3 + 0.5) + tgamma((l + 1) / 3 + 0.5) + tgamma((l + 2) / 3 + 0.5); // floor((l+i)/3) via int div
        return sqrt(numr / denm);
    };
    const double I_ab = 2 * M_PI * B_Ylm(a, r) * B_Ylm(b, s);
    double result = 0.0;
    for(int t = 0; t <= r / 2; t++){
        for(int u = 0; u <= s / 2; u++){
            result += E_rstu(r, s, t, u, a, b, QAB) * fact(t + u) / intpow(a + b, 1 + t + u);
        }
    }
    return I_ab * result;
}

double E_rstu(int r, int s, int t, int u, double a, double b, double QAB){
    return nCk(r/2, t) * nCk(s/2, u) * Theta(r/2 - t, s/2 - u, a, b, QAB);
}

double Theta(int k, int l, double a, double b, double QAB){
    auto d_exp2_d_expi = [](double exp, int i){
        assert(i >= 0);
        return (i != 0 ? (i != 1 ? (i != 2 ? (0.0) : 2.0) : 2.0 * exp) : exp * exp);
    };
    if((k < 0) || (l < 0)){return 0.0;}
    else if((k==0) && (l==0)){return exp(-a * b * QAB * QAB / (a + b));}
    else if(k==0){
        double sum = 0.0;
        for(int q = 0; q <= l - 1; q++){
            sum += nCk(l - 1, q) * Theta(0, l - 1 - q, a, b, QAB) * fact(q + 1) * a * a / intpow(a + b, q + 2);
        }
        return QAB * QAB * sum;
    }
    else{
        double sum = 0.0; 
        double dterm;
        for(int p = 0; p <= k - 1; p++){ 
            for(int q = 0; q <= l; q++){
                dterm = 0.0;
                for(int i = 0; i <= (q > 2 ? 2 : q); i++){
                    dterm += nCk(q, i) * ((i % 2 == 0) ? 1 : -1) * fact(p + q - i + 1) * d_exp2_d_expi(b, i) / intpow(a + b, 2 - i);
                }
                dterm /= intpow(a + b, p + q);
                sum += nCk(k - 1, p) * nCk(l, q) * Theta(k - 1 - p, l - q, a, b, QAB) * dterm;
            }
        }
        return QAB * QAB * sum;
    }
}
