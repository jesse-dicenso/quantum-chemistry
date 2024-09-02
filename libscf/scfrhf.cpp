#include "scfrhf.hpp"

Matrix R_density_matrix(Matrix C, int N){
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

Matrix R_F(Matrix Hcore, Matrix P, std::vector<std::vector<std::vector<std::vector<double>>>> g){
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

double R_E0(Matrix P, Matrix Hcore, Matrix F){
	double sum = 0;
	for(int i = 0; i < P.rows; i++){
		for(int j = 0; j < P.cols; j++){
			sum += P.matrix[j][i]*(Hcore.matrix[i][j]+F.matrix[i][j]);
		}
	}
	return 0.5 * sum;
}

void R_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N){
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

void R_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps, int* icd){
	// Uses commutation of F and P for error metric
	// Perform fixed-point iterations until iteration number
	// equals the subspace size, then perform DIIS iterations
	if(i < sps){
		Matrix d = *p;
		R_FPI(s, hcore, eris, x, p, f, fo, e, co, c, Eo, err, N);
		for(int j = 0; j < sps-1; j++){
			SPf[j] = SPf[j+1];
			SPe[j] = SPe[j+1];
		}
		SPf[sps-1] = *f;
		Matrix ev = transpose(x) * ((*f) * (*p) * s - s * (*p) * (*f)) * x;
		SPe[sps-1] = ev;
	}
	else{
		std::vector<Matrix> tec(2);
		std::vector<double> weights(sps+1);
		double tEo;
		double errmax;
	
		*f   = R_F(hcore, *p, eris);
		
		// Store f, update error vectors
		for(int j = 0; j < sps-1; j++){
			SPf[j] = SPf[j+1];
			SPe[j] = SPe[j+1];
		}
		SPf[sps-1] = *f;
		Matrix ev = transpose(x) * ((*f) * (*p) * s - s * (*p) * (*f)) * x;
		SPe[sps-1] = ev;

		errmax = SPe[0].matrix[0][0];
		for(int j = 0; j < sps; j++){
			for(int k = 0; k < SPe[j].rows; k++){
				for(int l = 0; l < SPe[j].cols; l++){
					if(fabs(SPe[j].matrix[k][l]) > fabs(errmax)){
						errmax = SPe[j].matrix[k][l];
					}
				}
			}
		}

		// Set up linear system and solve for weights
		Matrix LHS(sps+1, sps+1);
		for(int j = 0; j < LHS.rows; j++){
			for(int k = 0; k < LHS.cols; k++){
				if((j==LHS.rows-1) && (k==LHS.cols-1)){
					LHS.matrix[j][k] = 0;
				}
				else if((j==LHS.rows-1) || (k==LHS.cols-1)){
					LHS.matrix[j][k] = -1;
				}
				else{
					LHS.matrix[j][k] = Tr(SPe[j]*SPe[k]);
				}
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
		
		*Eo  = R_E0(*p, hcore, *f);
		*err = errmax;
		
		// Build fock matrix from previous fock matrices and weights
		*f = zero(p->rows, p->cols);
		for(int j = 0; j < sps; j++){
			*f = *f + SPf[j] * weights[j];
		}

		*fo  = transpose(x) * (*f) * x;
		tec  = diagonalize(*fo);
		*e   = tec[0];
		*co  = tec[1];
		*c   = x * (*co);
		*p   = R_density_matrix(*c, N);
	}
}
