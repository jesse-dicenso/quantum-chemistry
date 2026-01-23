#include "scfr.hpp"

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

void R_FPI(const Matrix& s, const Matrix& hcore, const Matrix& x, XC* xc, Matrix* f, Matrix* fo, 
		   Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i)
{
	// Uses Ediff between iterations for error metric
	std::vector<Matrix> tec(2);
	Matrix j = coulomb(*(xc->P), *(xc->eris));
	call_xc_functional(xc);

	*f = fock(hcore, j, *(xc->FXC));
	double tEo  = E0(*xc, hcore, j);
	*err = tEo - *Eo;
	*Eo = tEo;

	*fo  = transpose(x) * (*f) * x;
	tec  = diagonalize(*fo);
	*e   = tec[0];
	*co  = tec[1];
	*c   = x * (*co);
	*(xc->P) = R_density_matrix(*c, N);
	
	std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err << std::setw(10) << "fp" << std::endl;
}

void R_DIIS(const Matrix& s, const Matrix& hcore, const Matrix& x, XC* xc, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, 
			double* Eo, double* err, int N, int i, std::vector<Matrix>& SPf, std::vector<Matrix>& SPe, int sps, int* icd)
{
	// Uses commutation of F and P for error metric
	if(i <= 3){	
		R_FPI(s, hcore, x, xc, f, fo, e, co, c, Eo, err, N, i);
		SPf.push_back(*f);
		SPe.push_back((*f) * (*(xc->P)) * s - s * (*(xc->P)) * (*f));
	}
	else if(i < sps){
		std::vector<Matrix> tec(2);
		std::vector<double> weights(sps+1);
		double tEo;
		double terr = 0;
	
		Matrix j = coulomb(*(xc->P), *(xc->eris));
		call_xc_functional(xc);
		*f = fock(hcore, j, *(xc->FXC));
		
		SPf.push_back(*f);
		SPe.push_back((*f) * (*(xc->P)) * s - s * (*(xc->P)) * (*f));

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
			std::cout << "*** WARNING: DIIS SYSTEM ILL-CONDITIONED, SWITCHING TO FPI (min 3 iter) ***" << std::endl;
			return;
		}
		
		// Build fock matrix from previous fock matrices and weights
		*f = zero(xc->P->rows, xc->P->cols);
		for(int j = 0; j < n; j++){
			*f = *f + SPf[j] * weights[j];
		}

		*Eo  = E0(*xc, hcore, j);

		*fo  = transpose(x) * (*f) * x;
		tec  = diagonalize(*fo);
		*e   = tec[0];
		*co  = tec[1];
		*c   = x * (*co);
		*(xc->P) = R_density_matrix(*c, N);	

		std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err;
		std::cout << std::setw(9) << "diis(" << SPe.size() << ")" << std::endl;
	}
	else{
		std::vector<Matrix> tec(2);
		std::vector<double> weights(sps+1);
		double tEo;
		double terr = 0;

		Matrix j = coulomb(*(xc->P), *(xc->eris));
		call_xc_functional(xc);
		*f = fock(hcore, j, *(xc->FXC));
		
		SPf.push_back(*f);
		SPe.push_back((*f) * (*(xc->P)) * s - s * (*(xc->P)) * (*f));

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
			std::cout << "*** WARNING: DIIS SYSTEM ILL-CONDITIONED, SWITCHING TO FPI (min 3 iter) ***" << std::endl;
			return;
		}
		
		// Build fock matrix from previous fock matrices and weights
		*f = zero(xc->P->rows, xc->P->cols);
		for(int j = 0; j < sps; j++){
			*f = *f + SPf[j] * weights[j];
		}
	
		*Eo  = E0(*xc, hcore, j);

		*fo  = transpose(x) * (*f) * x;
		tec  = diagonalize(*fo);
		*e   = tec[0];
		*co  = tec[1];
		*c   = x * (*co);
		*(xc->P) = R_density_matrix(*c, N);

		std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err;
		std::cout << std::setw(9) << "diis(" << sps << ")" << std::endl;

		SPf.erase(SPf.begin());		
		SPe.erase(SPe.begin());		
	}
}
