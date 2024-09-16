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

double UR_E0(Matrix PT, Matrix Pa, Matrix Pb, Matrix Hcore, Matrix Fa, Matrix Fb){
	double sum = 0;
	for(int i = 0; i < PT.rows; i++){
		for(int j = 0; j < PT.cols; j++){
			sum += PT.matrix[j][i] * Hcore.matrix[i][j] + Pa.matrix[j][i] * Fa.matrix[i][j] + Pb.matrix[j][i] * Fb.matrix[i][j];
		}
	}
	return 0.5 * sum;
}

void UR_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb){
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

void UR_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, Matrix* SPfa, Matrix* SPfb, Matrix* SPea, Matrix* SPeb, int sps, int* icd){
	// Uses commutation of F and P for error metric
	// Perform fixed-point iterations until iteration number
	// equals the subspace size, then perform DIIS iterations
	if(i < sps){
		Matrix da = *pa;
		Matrix db = *pb;
		UR_FPI(s, hcore, eris, x, pt, pa, pb, fa, fb, fao, fbo, ea, eb, cao, cbo, ca, cb, Eo, err, Na, Nb);
		for(int j = 0; j < sps-1; j++){
			SPfa[j] = SPfa[j+1];
			SPfb[j] = SPfb[j+1];
			SPea[j] = SPea[j+1];
			SPeb[j] = SPeb[j+1];
		}
		SPfa[sps-1] = *fa;
		SPfb[sps-1] = *fb;
		Matrix eva = transpose(x) * ((*fa) * (*pa) * s - s * (*pa) * (*fa)) * x;
		Matrix evb = transpose(x) * ((*fb) * (*pb) * s - s * (*pb) * (*fb)) * x;
		SPea[sps-1] = eva;
		SPeb[sps-1] = evb;
	}
	else{
		std::vector<Matrix> tec_a(2);
		std::vector<Matrix> tec_b(2);
		std::vector<double> weights_a(sps+1);
		std::vector<double> weights_b(sps+1);
		double errmax_a;
		double errmax_b;
	
		*fa = UR_F(hcore, *pt, *pa, eris);
		*fb = UR_F(hcore, *pt, *pb, eris);
		
		// Store fa, fb, update error vectors
		for(int j = 0; j < sps-1; j++){
			SPfa[j] = SPfa[j+1];
			SPfb[j] = SPfb[j+1];
			SPea[j] = SPea[j+1];
			SPeb[j] = SPeb[j+1];
		}
		SPfa[sps-1] = *fa;
		SPfb[sps-1] = *fb;
		Matrix eva = transpose(x) * ((*fa) * (*pa) * s - s * (*pa) * (*fa)) * x;
		Matrix evb = transpose(x) * ((*fb) * (*pb) * s - s * (*pb) * (*fb)) * x;
		SPea[sps-1] = eva;
		SPeb[sps-1] = evb;
		
		errmax_a = SPea[0].matrix[0][0];
		errmax_b = SPeb[0].matrix[0][0];

		for(int j = 0; j < sps; j++){
			for(int k = 0; k < SPea[j].rows; k++){
				for(int l = 0; l < SPea[j].cols; l++){
					if(fabs(SPea[j].matrix[k][l]) > fabs(errmax_a)){
						errmax_a = SPea[j].matrix[k][l];
					}
				}
			}
			for(int k = 0; k < SPeb[j].rows; k++){
				for(int l = 0; l < SPeb[j].cols; l++){
					if(fabs(SPeb[j].matrix[k][l]) > fabs(errmax_b)){
						errmax_b = SPeb[j].matrix[k][l];
					}
				}
			}
		}
		
		// Set up linear system and solve for weights
		Matrix LHS_a(sps+1, sps+1);
		for(int j = 0; j < LHS_a.rows; j++){
			for(int k = 0; k < LHS_a.cols; k++){
				if((j==LHS_a.rows-1) && (k==LHS_a.cols-1)){
					LHS_a.matrix[j][k] = 0;
				}
				else if((j==LHS_a.rows-1) || (k==LHS_a.cols-1)){
					LHS_a.matrix[j][k] = -1;
				}
				else{
					LHS_a.matrix[j][k] = Tr(SPea[j]*SPea[k]);
				}
			}
		}
		Matrix LHS_b(sps+1, sps+1);
		for(int j = 0; j < LHS_b.rows; j++){
			for(int k = 0; k < LHS_b.cols; k++){
				if((j==LHS_b.rows-1) && (k==LHS_b.cols-1)){
					LHS_b.matrix[j][k] = 0;
				}
				else if((j==LHS_b.rows-1) || (k==LHS_b.cols-1)){
					LHS_b.matrix[j][k] = -1;
				}
				else{
					LHS_b.matrix[j][k] = Tr(SPeb[j]*SPeb[k]);
				}
			}
		}

		Matrix RHS_a(sps+1, 1);
		RHS_a.matrix[sps][0] = -1;
		Matrix RHS_b(sps+1, 1);
		RHS_b.matrix[sps][0] = -1;

		weights_a = sym_linear_solve(LHS_a, RHS_a, icd);
		weights_b = sym_linear_solve(LHS_b, RHS_b, icd);
		
		if(*icd!=0){
			std::cout << "*** WARNING: DIIS SYSTEM ILL-CONDITIONED, SWITCHING TO FPI (min 3 iter) ***\n";
			std::cout.flush();
			return;
		}
		
		//double tEo = UR_E0(*pt, *pa, *pb, hcore, *fa, *fb);
		//*err = tEo - *Eo;
		*Eo = UR_E0(*pt, *pa, *pb, hcore, *fa, *fb);
		
		if(fabs(errmax_a) > fabs(errmax_b)){
			*err = errmax_a;
		}
		else{
			*err = errmax_b;
		}	
		
		// Build fock matrices from previous fock matrices and weights
		*fa = zero(pt->rows, pt->cols);
		*fb = zero(pt->rows, pt->cols);
		for(int j = 0; j < sps; j++){
			*fa = *fa + SPfa[j] * weights_a[j];
			*fb = *fb + SPfb[j] * weights_b[j];
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
}
