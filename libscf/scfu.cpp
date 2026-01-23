#include "scfu.hpp"

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

void UR_FPI(const Matrix& s, const Matrix& hcore, const Matrix& x, XC* xc, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, 
		    Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i)
{
	std::vector<Matrix> tec_a(2);
	std::vector<Matrix> tec_b(2);

	Matrix j  = coulomb(*(xc->P), *(xc->eris));
	call_xc_functional(xc);
	*fa = fock(hcore, j, *(xc->FXC_A));
	*fb = fock(hcore, j, *(xc->FXC_B));
	
	double tEo = E0(*xc, hcore, j);
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

	*(xc->P_A) = UR_density_matrix(*ca, Na);
	*(xc->P_B) = UR_density_matrix(*cb, Nb);
	*(xc->P)   = (*(xc->P_A)) + (*(xc->P_B));
	
	std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err << std::setw(10) << "fp" << std::endl;
}

void UR_DIIS(const Matrix& s, const Matrix& hcore, const Matrix& x, XC* xc, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, 
			 Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, 
			 std::vector<Matrix>& SPfa, std::vector<Matrix>& SPfb, std::vector<Matrix>& SPea, std::vector<Matrix>& SPeb, int sps, int* icd)
{
	// Uses commutation of F and P for error metric
	if(i <= 3){
		UR_FPI(s, hcore, x, xc, fa, fb, fao, fbo, ea, eb, cao, cbo, ca, cb, Eo, err, Na, Nb, i);
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*(xc->P_A)) * s - s * (*(xc->P_A)) * (*fa));
		SPeb.push_back((*fb) * (*(xc->P_B)) * s - s * (*(xc->P_B)) * (*fb));
	}
	else if(i < sps){
		std::vector<Matrix> tec_a(2);
		std::vector<Matrix> tec_b(2);
		std::vector<double> weights(sps+1);
		double terr = 0;
	
		Matrix j    = coulomb(*(xc->P), *(xc->eris));
		call_xc_functional(xc);
		*fa = fock(hcore, j, *(xc->FXC_A));
		*fb = fock(hcore, j, *(xc->FXC_B));
	
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*(xc->P_A)) * s - s * (*(xc->P_A)) * (*fa));
		SPeb.push_back((*fb) * (*(xc->P_B)) * s - s * (*(xc->P_B)) * (*fb));

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
			std::cout << "*** WARNING: DIIS SYSTEM ILL-CONDITIONED, SWITCHING TO FPI (min 3 iter) ***" << std::endl;
			return;
		}
		
		*Eo = E0(*xc, hcore, j);
		
		// Build fock matrices from previous fock matrices and weights
		*fa = zero(xc->P->rows, xc->P->cols);
		*fb = zero(xc->P->rows, xc->P->cols);
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
		*(xc->P_A) = UR_density_matrix(*ca, Na);
		*(xc->P_B) = UR_density_matrix(*cb, Nb);
		*(xc->P)   = *(xc->P_A) + *(xc->P_B);
		
		std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err;
		std::cout << std::setw(9) << "diis(" << SPea.size() << ")" << std::endl;
	}
	else{
		std::vector<Matrix> tec_a(2);
		std::vector<Matrix> tec_b(2);
		std::vector<double> weights(sps+1);
		double terr = 0;
	
		Matrix j  = coulomb(*(xc->P), *(xc->eris));
		call_xc_functional(xc);
		*fa = fock(hcore, j, *(xc->FXC_A));
		*fb = fock(hcore, j, *(xc->FXC_B));
	
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*(xc->P_A)) * s - s * (*(xc->P_A)) * (*fa));
		SPeb.push_back((*fb) * (*(xc->P_B)) * s - s * (*(xc->P_B)) * (*fb));
	
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
			std::cout << "*** WARNING: DIIS SYSTEM ILL-CONDITIONED, SWITCHING TO FPI (min 3 iter) ***" << std::endl;
			return;
		}
		
		*Eo = E0(*xc, hcore, j);
		
		// Build fock matrices from previous fock matrices and weights
		*fa = zero(xc->P->rows, xc->P->cols);
		*fb = zero(xc->P->rows, xc->P->cols);
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
		*(xc->P_A) = UR_density_matrix(*ca, Na);
		*(xc->P_B) = UR_density_matrix(*cb, Nb);
		*(xc->P)   = *(xc->P_A) + *(xc->P_B);
		
		std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err;
		std::cout << std::setw(9) << "diis(" << sps << ")" << std::endl;

		SPfa.erase(SPfa.begin());
		SPfb.erase(SPfb.begin());
		SPea.erase(SPea.begin());
		SPeb.erase(SPeb.begin());
	}
}
