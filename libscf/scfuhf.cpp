#include "scfuhf.hpp"

Matrix UR_density_matrix(const Matrix& C, int N){
	Matrix P(C.rows, C.cols);
	for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			double sum = 0;
			for(int a = 0; a < N; a++){
				sum += C.matrix[i][a]*C.matrix[j][a];
			} 
			P.matrix[i][j] = sum;
		}
	}
	return P;
}

Matrix UR_F(const Matrix& Hcore, const Matrix& PT, const Matrix& Ps, const std::vector<std::vector<std::vector<std::vector<double>>>>& g){
	Matrix F(PT.rows, PT.cols);
	double sum;
	for(int mu = 0; mu < F.rows; mu++){
		for(int nu = 0; nu < F.rows; nu++){
			sum = Hcore.matrix[mu][nu];
			for(int ld = 0; ld < F.rows; ld++){
				for(int sg = 0; sg < F.rows; sg++){
					sum += (PT.matrix[ld][sg] * g[mu][nu][sg][ld]) - (Ps.matrix[ld][sg] * g[mu][ld][sg][nu]);
				}
			}
			F.matrix[mu][nu] = sum;
		}
	}
	return F;
}

double UR_E0(const Matrix& PT, const Matrix& Pa, const Matrix& Pb, const Matrix& Hcore, const Matrix& Fa, const Matrix& Fb){
	double sum = 0;
	for(int i = 0; i < PT.rows; i++){
		for(int j = 0; j < PT.cols; j++){
			sum += PT.matrix[j][i] * Hcore.matrix[i][j] + Pa.matrix[j][i] * Fa.matrix[i][j] + Pb.matrix[j][i] * Fb.matrix[i][j];
		}
	}
	return 0.5 * sum;
}

void UR_FPI(const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, const Matrix& x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb){
	std::vector<Matrix> tec_a(2);
	std::vector<Matrix> tec_b(2);
	double tEo;

	*fa = UR_F(hcore, *pt, *pa, eris);
	*fb = UR_F(hcore, *pt, *pb, eris);
	
	tEo = UR_E0(*pt, *pa, *pb, hcore, *fa, *fb);
	*err = tEo - *Eo;
	*Eo = tEo;

	*fao  = transpose(x) * (*fa) * x;
	*fbo  = transpose(x) * (*fb) * x;
	tec_a  = diagonalize(*fao);
	tec_b  = diagonalize(*fbo);
	*ea  = tec_a[0];
	*cao = tec_a[1];
	*eb  = tec_b[0];
	*cbo = tec_b[1];
	*ca  = x * (*cao);
	*cb  = x * (*cbo);

	*pa  = UR_density_matrix(*ca, Na);
	*pb  = UR_density_matrix(*cb, Nb);
	*pt  = (*pa) + (*pb);
}

void UR_DIIS(const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, const Matrix& x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, std::vector<Matrix>& SPfa, std::vector<Matrix>& SPfb, std::vector<Matrix>& SPea, std::vector<Matrix>& SPeb, int sps, int* icd){
	// Uses commutation of F and P for error metric
	if(i <= 3){
		UR_FPI(s, hcore, eris, x, pt, pa, pb, fa, fb, fao, fbo, ea, eb, cao, cbo, ca, cb, Eo, err, Na, Nb);
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*pa) * s - s * (*pa) * (*fa));
		SPeb.push_back((*fb) * (*pb) * s - s * (*pb) * (*fb));
	}
	else if(i < sps){
		std::vector<Matrix> tec_a(2);
		std::vector<Matrix> tec_b(2);
		std::vector<double> weights(sps+1);
		double terr = 0;
	
		*fa = UR_F(hcore, *pt, *pa, eris);
		*fb = UR_F(hcore, *pt, *pb, eris);
		
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*pa) * s - s * (*pa) * (*fa));
		SPeb.push_back((*fb) * (*pb) * s - s * (*pb) * (*fb));

		int n = SPea.size();
	
		for(int j = 0; j < SPea.back().rows; j++){
			for(int k = 0; k < SPea.back().cols; k++){
				terr += SPea.back().matrix[j][k] * SPea.back().matrix[j][k] + 
						SPeb.back().matrix[j][k] * SPeb.back().matrix[j][k];
			}
		}

		*err = sqrt(terr);
	
		// Set up linear system and solve for weights
		Matrix LHS(n+1, n+1);
		LHS.matrix[n][n] = 0;
		for(int j = 0; j < n; j++){
			LHS.matrix[n][j] = -1;
			LHS.matrix[j][n] = -1;
		}
		for(int j = 0; j < n; j++){
			for(int k = 0; k < n; k++){
				LHS.matrix[j][k] = Tr(SPea[j] * SPea[k] + SPeb[j] * SPeb[k]);
			}
		}

		Matrix RHS(n+1, 1);
		RHS.matrix[n][0] = -1;

		weights = sym_linear_solve(LHS, RHS, icd);
		
		if(*icd!=0){
			std::cout << "*** WARNING: DIIS SYSTEM ILL-CONDITIONED, SWITCHING TO FPI (min 3 iter) ***\n";
			std::cout.flush();
			return;
		}
		
		*Eo = UR_E0(*pt, *pa, *pb, hcore, *fa, *fb);
		
		// Build fock matrices from previous fock matrices and weights
		*fa = zero(pt->rows, pt->cols);
		*fb = zero(pt->rows, pt->cols);
		for(int j = 0; j < n; j++){
			*fa = *fa + SPfa[j] * weights[j];
			*fb = *fb + SPfb[j] * weights[j];
		}

		*fao  = transpose(x) * (*fa) * x;
		*fbo  = transpose(x) * (*fb) * x;
		tec_a = diagonalize(*fao);
		tec_b = diagonalize(*fbo);
		*ea   = tec_a[0];
		*eb   = tec_b[0];
		*cao  = tec_a[1];
		*cbo  = tec_b[1];
		*ca   = x * (*cao);
		*cb   = x * (*cbo);
		*pa   = UR_density_matrix(*ca, Na);
		*pb   = UR_density_matrix(*cb, Nb);
		*pt   = *pa + *pb;
	}
	else{
		std::vector<Matrix> tec_a(2);
		std::vector<Matrix> tec_b(2);
		std::vector<double> weights(sps+1);
		double terr = 0;
	
		*fa = UR_F(hcore, *pt, *pa, eris);
		*fb = UR_F(hcore, *pt, *pb, eris);
		
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*pa) * s - s * (*pa) * (*fa));
		SPeb.push_back((*fb) * (*pb) * s - s * (*pb) * (*fb));
	
		for(int j = 0; j < SPea.back().rows; j++){
			for(int k = 0; k < SPea.back().cols; k++){
				terr += SPea.back().matrix[j][k] * SPea.back().matrix[j][k] + 
						SPeb.back().matrix[j][k] * SPeb.back().matrix[j][k];
			}
		}

		*err = sqrt(terr);
	
		// Set up linear system and solve for weights
		Matrix LHS(sps+1, sps+1);
		LHS.matrix[sps][sps] = 0;
		for(int j = 0; j < sps; j++){
			LHS.matrix[sps][j] = -1;
			LHS.matrix[j][sps] = -1;
		}
		for(int j = 0; j < sps; j++){
			for(int k = 0; k < sps; k++){
				LHS.matrix[j][k] = Tr(SPea[j] * SPea[k] + SPeb[j] * SPeb[k]);
			}
		}

		Matrix RHS(sps+1, 1);
		RHS.matrix[sps][0] = -1;

		weights = sym_linear_solve(LHS, RHS, icd);
		
		if(*icd!=0){
			std::cout << "*** WARNING: DIIS SYSTEM ILL-CONDITIONED, SWITCHING TO FPI (min 3 iter) ***\n";
			std::cout.flush();
			return;
		}
		
		*Eo = UR_E0(*pt, *pa, *pb, hcore, *fa, *fb);
		
		// Build fock matrices from previous fock matrices and weights
		*fa = zero(pt->rows, pt->cols);
		*fb = zero(pt->rows, pt->cols);
		for(int j = 0; j < sps; j++){
			*fa = *fa + SPfa[j] * weights[j];
			*fb = *fb + SPfb[j] * weights[j];
		}

		*fao  = transpose(x) * (*fa) * x;
		*fbo  = transpose(x) * (*fb) * x;
		tec_a = diagonalize(*fao);
		tec_b = diagonalize(*fbo);
		*ea   = tec_a[0];
		*eb   = tec_b[0];
		*cao  = tec_a[1];
		*cbo  = tec_b[1];
		*ca   = x * (*cao);
		*cb   = x * (*cbo);
		*pa   = UR_density_matrix(*ca, Na);
		*pb   = UR_density_matrix(*cb, Nb);
		*pt   = *pa + *pb;

		SPfa.erase(SPfa.begin());
		SPfb.erase(SPfb.begin());
		SPea.erase(SPea.begin());
		SPeb.erase(SPeb.begin());
	}
}
