#include "scfrhf.hpp"

Matrix R_density_matrix(const Matrix& C, int N){
	Matrix P(C.rows, C.cols);
	for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			double sum = 0;
			for(int a = 0; a < N/2; a++){
				sum += C.matrix[i][a]*C.matrix[j][a];
			} 
			P.matrix[i][j] = 2*sum;
		}
	}
	return P;
}

Matrix R_F(const Matrix& Hcore,const Matrix& P, const std::vector<std::vector<std::vector<std::vector<double>>>>& g){
	Matrix G(P.rows, P.cols);
	double sum;
	for(int mu = 0; mu < G.rows; mu++){
		for(int nu = 0; nu < G.rows; nu++){
			sum = 0;
			for(int ld = 0; ld < G.rows; ld++){
				for(int sg = 0; sg < G.rows; sg++){
					sum += P.matrix[ld][sg] * (g[mu][nu][sg][ld] - 0.5*g[mu][ld][sg][nu]);
				}
			}
			G.matrix[mu][nu] = sum;
		}
	}
	return Hcore + G;
}

double R_E0(const Matrix& P, const Matrix& Hcore, const Matrix& F){
	double sum = 0;
	for(int i = 0; i < P.rows; i++){
		for(int j = 0; j < P.cols; j++){
			sum += P.matrix[j][i]*(Hcore.matrix[i][j]+F.matrix[i][j]);
		}
	}
	return 0.5 * sum;
}

void R_FPI(const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, const Matrix& x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N){
	// Uses Ediff between iterations for error metric
	std::vector<Matrix> tec(2);
	double tEo;

	*f = R_F(hcore, *p, eris);
	tEo  = R_E0(*p, hcore, *f);
	*err = tEo - *Eo;
	*Eo = tEo;

	*fo  = transpose(x) * (*f) * x;
	tec  = diagonalize(*fo);
	*e   = tec[0];
	*co  = tec[1];
	*c   = x * (*co);
	*p   = R_density_matrix(*c, N);
}

void R_DIIS(const Matrix& s, const Matrix& hcore, const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, const Matrix& x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, std::vector<Matrix>& SPf, std::vector<Matrix>& SPe, int sps, int* icd){
	// Uses commutation of F and P for error metric
	if(i == 1){	
		R_FPI(s, hcore, eris, x, p, f, fo, e, co, c, Eo, err, N);
		SPf.push_back(*f);
		Matrix ev = ((*f) * (*p) * s - s * (*p) * (*f));
		SPe.push_back(ev);	
	}
	else if(i < sps){
		std::vector<Matrix> tec(2);
		std::vector<double> weights(sps+1);
		double tEo;
		double terr = 0;
	
		*f   = R_F(hcore, *p, eris);
		
		SPf.push_back(*f);
		Matrix ev = ((*f) * (*p) * s - s * (*p) * (*f));
		SPe.push_back(ev);

		int n = SPe.size();

		for(int j = 0; j < SPe.back().rows; j++){
			for(int k = 0; k < SPe.back().cols; k++){
				terr += SPe.back().matrix[j][k] * SPe.back().matrix[j][k];
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
				LHS.matrix[j][k] = Tr(SPe[j]*SPe[k]);
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
		
		// Build fock matrix from previous fock matrices and weights
		*f = zero(p->rows, p->cols);
		for(int j = 0; j < n; j++){
			*f = *f + SPf[j] * weights[j];
		}

		*Eo  = R_E0(*p, hcore, *f);

		*fo  = transpose(x) * (*f) * x;
		tec  = diagonalize(*fo);
		*e   = tec[0];
		*co  = tec[1];
		*c   = x * (*co);
		*p   = R_density_matrix(*c, N);	
	}
	else{
		std::vector<Matrix> tec(2);
		std::vector<double> weights(sps+1);
		double tEo;
		double terr = 0;

		*f   = R_F(hcore, *p, eris);
		
		SPf.push_back(*f);
		Matrix ev = ((*f) * (*p) * s - s * (*p) * (*f));
		SPe.push_back(ev);

		for(int j = 0; j < SPe.back().rows; j++){
			for(int k = 0; k < SPe.back().cols; k++){
				terr += SPe.back().matrix[j][k] * SPe.back().matrix[j][k];
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
				LHS.matrix[j][k] = Tr(SPe[j]*SPe[k]);
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
		
		// Build fock matrix from previous fock matrices and weights
		*f = zero(p->rows, p->cols);
		for(int j = 0; j < sps; j++){
			*f = *f + SPf[j] * weights[j];
		}
		
		*Eo  = R_E0(*p, hcore, *f);

		*fo  = transpose(x) * (*f) * x;
		tec  = diagonalize(*fo);
		*e   = tec[0];
		*co  = tec[1];
		*c   = x * (*co);
		*p   = R_density_matrix(*c, N);

		SPf.erase(SPf.begin());		
		SPe.erase(SPe.begin());		
	}
}
