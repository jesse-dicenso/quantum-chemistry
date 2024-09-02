#include "scfuhf.hpp"

Matrix UR_density_matrix(Matrix C, int N){
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

Matrix UR_F(Matrix Hcore, Matrix PT, Matrix Ps, std::vector<std::vector<std::vector<std::vector<double>>>> g){
	Matrix G(PT.rows, PT.cols);
	double sum;
	for(int mu = 0; mu < G.rows; mu++){
		for(int nu = 0; nu < G.rows; nu++){
			sum = 0;
			for(int ld = 0; ld < G.rows; ld++){
				for(int sg = 0; sg < G.rows; sg++){
					sum += PT.matrix[ld][sg] * g[mu][nu][sg][ld] - Ps.matrix[ld][sg] * g[mu][ld][sg][nu];
				}
			}
			G.matrix[mu][nu] = sum;
		}
	}
	return Hcore + G;
}

double UR_E0(Matrix PT, Matrix Pa, Matrix Pb, Matrix Hcore, Matrix Fa, Matrix Fb){
	double sum = 0;
	for(int i = 0; i < PT.rows; i++){
		for(int j = 0; j < PT.cols; j++){
			sum += PT.matrix[j][i] * Hcore.matrix[i][j] + Pa.matrix[j][i] * Fa.matrix[i][j] + Pb.matrix[j][i] * Fb.matrix[i][j];
		}
	}
	return 0.5 * sum;
}

void UR_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int N){
	std::vector<Matrix> tec_a(2);
	std::vector<Matrix> tec_b(2);
	double tEo;

	*fa = UR_F(hcore, *pt, *pb, eris);
	*fb = UR_F(hcore, *pt, *pa, eris);
	
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
	*pa  = UR_density_matrix(*ca, N);
	*pb  = UR_density_matrix(*cb, N);
	*pt  = (*pa) + (*pb);
}
/*
void UR_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* foa, Matrix* fob, Matrix* ea, Matrix* eb, Matrix* coa, Matrix* cob, Matrix* ca, Matrix* cb, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps){
	// Perform fixed-point iterations until iteration number
	// equals the subspace size, then perform DIIS iterations
	if(i < sps){
		Matrix d = *p;
		FPI(R, s, hcore, eris, x, p, g, f, fo, e, co, c, Eo, err, N);
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
	
		*g   = G(*p, eris);
		*f   = hcore + *g;
		*Eo  = E0(*p, hcore, *f);
		
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
		*err = errmax;

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
		weights = sym_linear_solve(LHS, RHS);
		
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
		*p   = UR_density_matrix(*c, N);
	}
}
*/
