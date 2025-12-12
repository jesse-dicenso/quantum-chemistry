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

double UR_E0(const Matrix& PT, const Matrix& Pa, const Matrix& Pb, const Matrix& Hcore, const Matrix& Fa, const Matrix& Fb){
	double sum = 0;
	for(int i = 0; i < PT.rows; i++){
		for(int j = 0; j < PT.cols; j++){
			sum += PT.matrix[j][i] * Hcore.matrix[i][j] + Pa.matrix[j][i] * Fa.matrix[i][j] + Pb.matrix[j][i] * Fb.matrix[i][j];
		}
	}
	return 0.5 * sum;
}

void UR_FPI (const Matrix& s, const Matrix& hcore, const Matrix& x, XC_inp* xc_inp, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, 
			 Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i)
{
	std::vector<Matrix> tec_a(2);
	std::vector<Matrix> tec_b(2);
	double tEo;

	Matrix j  = coulomb(*(xc_inp->PT), *(xc_inp->eris));
	XC_ret fxc = F_XC(xc_inp);
	Matrix fxca = fxc.F_XC_1;
	Matrix fxcb = fxc.F_XC_2;

	*fa = fock(hcore, j, fxca);
	*fb = fock(hcore, j, fxcb);
	
	tEo = UR_E0(*(xc_inp->PT), *(xc_inp->PA), *(xc_inp->PB), hcore, *fa, *fb);
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

	*(xc_inp->PA) = UR_density_matrix(*ca, Na);
	*(xc_inp->PB) = UR_density_matrix(*cb, Nb);
	*(xc_inp->PT) = (*(xc_inp->PA)) + (*(xc_inp->PB));
	
	std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err << std::setw(10) << "fp" << std::endl;
}

void UR_DIIS(const Matrix& s, const Matrix& hcore, const Matrix& x, XC_inp* xc_inp, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, 
			 Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, 
			 std::vector<Matrix>& SPfa, std::vector<Matrix>& SPfb, std::vector<Matrix>& SPea, std::vector<Matrix>& SPeb, int sps, int* icd)
{
	// Uses commutation of F and P for error metric
	if(i <= 3){
		UR_FPI(s, hcore, x, xc_inp, fa, fb, fao, fbo, ea, eb, cao, cbo, ca, cb, Eo, err, Na, Nb, i);
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*(xc_inp->PA)) * s - s * (*(xc_inp->PA)) * (*fa));
		SPeb.push_back((*fb) * (*(xc_inp->PB)) * s - s * (*(xc_inp->PB)) * (*fb));
	}
	else if(i < sps){
		std::vector<Matrix> tec_a(2);
		std::vector<Matrix> tec_b(2);
		std::vector<double> weights(sps+1);
		double terr = 0;
	
		Matrix j    = coulomb(*(xc_inp->PT), *(xc_inp->eris));
		XC_ret fxc  = F_XC(xc_inp);
		Matrix fxca = fxc.F_XC_1;
		Matrix fxcb = fxc.F_XC_2;
	
		*fa = fock(hcore, j, fxca);
		*fb = fock(hcore, j, fxcb);
	
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*(xc_inp->PA)) * s - s * (*(xc_inp->PA)) * (*fa));
		SPeb.push_back((*fb) * (*(xc_inp->PB)) * s - s * (*(xc_inp->PB)) * (*fb));

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
		
		*Eo = UR_E0(*(xc_inp->PT), *(xc_inp->PA), *(xc_inp->PB), hcore, *fa, *fb);
		
		// Build fock matrices from previous fock matrices and weights
		*fa = zero(xc_inp->PT->rows, xc_inp->PT->cols);
		*fb = zero(xc_inp->PT->rows, xc_inp->PT->cols);
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
		*(xc_inp->PA) = UR_density_matrix(*ca, Na);
		*(xc_inp->PB) = UR_density_matrix(*cb, Nb);
		*(xc_inp->PT) = *(xc_inp->PA) + *(xc_inp->PB);
		
		std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err;
		std::cout << std::setw(9) << "diis(" << SPea.size() << ")" << std::endl;
	}
	else{
		std::vector<Matrix> tec_a(2);
		std::vector<Matrix> tec_b(2);
		std::vector<double> weights(sps+1);
		double terr = 0;
	
		Matrix j  = coulomb(*(xc_inp->PT), *(xc_inp->eris));
		XC_ret fxc = F_XC(xc_inp);
		Matrix fxca = fxc.F_XC_1;
		Matrix fxcb = fxc.F_XC_2;
	
		*fa = fock(hcore, j, fxca);
		*fb = fock(hcore, j, fxcb);
	
		SPfa.push_back(*fa);
		SPfb.push_back(*fb);
		SPea.push_back((*fa) * (*(xc_inp->PA)) * s - s * (*(xc_inp->PA)) * (*fa));
		SPeb.push_back((*fb) * (*(xc_inp->PB)) * s - s * (*(xc_inp->PB)) * (*fb));
	
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
		
		*Eo = UR_E0(*(xc_inp->PT), *(xc_inp->PA), *(xc_inp->PB), hcore, *fa, *fb);
		
		// Build fock matrices from previous fock matrices and weights
		*fa = zero(xc_inp->PT->rows, xc_inp->PT->cols);
		*fb = zero(xc_inp->PT->rows, xc_inp->PT->cols);
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
		*(xc_inp->PA) = UR_density_matrix(*ca, Na);
		*(xc_inp->PB) = UR_density_matrix(*cb, Nb);
		*(xc_inp->PT) = *(xc_inp->PA) + *(xc_inp->PB);
		
		std::cout << std::setw(3) << i << std::setw(20) << *Eo << std::setw(20) << *err;
		std::cout << std::setw(9) << "diis(" << sps << ")" << std::endl;

		SPfa.erase(SPfa.begin());
		SPfb.erase(SPfb.begin());
		SPea.erase(SPea.begin());
		SPeb.erase(SPeb.begin());
	}
}
