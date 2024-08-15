#include "linalg.hpp"

Matrix::Matrix(int r, int c, bool sym){
	rows = r;
	cols = c;
	isSymmetric = sym;
	matrix = new double *[rows];
	for(int i = 0; i < rows; i++){
		matrix[i] = new double[cols];
		for(int j = 0; j < cols; j++){
			matrix[i][j] = 0;
		}
	}
}

Matrix::Matrix(const Matrix& A){
	rows = A.rows;
	cols = A.cols;
	isSymmetric = A.isSymmetric;
	matrix = new double *[rows];
	for(int i = 0; i < rows; i++){
		matrix[i] = new double[cols];
		for(int j = 0; j < cols; j++){
			matrix[i][j] = A.matrix[i][j];
		}
	}
}

Matrix::Matrix(){
	rows = 1;
	cols = 1;
	matrix = new double *[rows];
	matrix[0] = new double[cols];
}

Matrix::~Matrix(){
	for(int i = 0; i < rows; i++){
		delete[] matrix[i];
	}
	delete[] matrix;
}

void Matrix::printMatrix(){
	for(int i = 0; i < rows; i++){
		std::cout << "     ";
		for(int j = 0; j < cols; j++){
			std::cout << matrix[i][j] << std::setw(20);
		}
		std::cout << '\n';
	}
}

Matrix Matrix::operator-() const{
	Matrix mat(rows, cols);
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			mat.matrix[i][j] == -matrix[i][j];
		}
	}
	return mat;
}

Matrix& Matrix::operator=(const Matrix& A){
	rows = A.rows;
	cols = A.cols;
	isSymmetric = A.isSymmetric;	
	for(int h = 0; h < rows; h++){
		delete[] matrix[h];
	}
	delete[] matrix;
	
	matrix = new double *[rows];
	for(int i = 0; i < rows; i++){
		matrix[i] = new double[cols];
		for(int j = 0; j < cols; j++){
			matrix[i][j] = A.matrix[i][j];
		}
	}
	return *this;
}

Matrix Matrix::operator+(const Matrix& A) const{
	assert((rows==A.rows) && (cols==A.cols));
	Matrix sum(rows, cols);
	if(A.isSymmetric){
		sum.isSymmetric=true;
	}
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			sum.matrix[i][j] = matrix[i][j] + A.matrix[i][j];
		}
	}
	return sum;
}

Matrix Matrix::operator-(const Matrix& A) const{
	assert((rows==A.rows) && (cols==A.cols));
	Matrix dif(rows, cols);
	if(A.isSymmetric){
		dif.isSymmetric=true;
	}
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			dif.matrix[i][j] = matrix[i][j] - A.matrix[i][j];
		}
	}
	return dif;
}

Matrix Matrix::operator*(double c) const{
	Matrix mat(rows, cols);
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			mat.matrix[i][j] = matrix[i][j] * c;
		}	
	}
	return mat;
}

Matrix Matrix::operator*(const Matrix& A) const{
	assert(cols==A.rows);
	Matrix product(rows, A.cols);
	for(int i = 0; i < product.rows; i++){
		for(int j = 0; j < product.cols; j++){
			for(int k = 0; k < cols; k++){
				product.matrix[i][j] += matrix[i][k] * A.matrix[k][j];
			}
		}
	}
	return product;	
}

Matrix transpose(const Matrix A){
	Matrix tp(A.cols, A.rows);
	for(int i = 0; i < A.rows; i++){
		for(int j = 0; j < A.cols; j++){
			tp.matrix[j][i] = A.matrix[i][j];
		}
	}
	return tp;
}

Matrix Matrix::getrow(int i){
	Matrix M(1, cols);
	for(int j = 0; j < cols; j++){
		M.matrix[0][j] = matrix[i][j];
	}
	return M;
}

Matrix Matrix::getcol(int i){
	Matrix M(rows, 1);
	for(int j = 0; j < rows; j++){
		M.matrix[j][0] = matrix[j][i];
	}
	return M;
}

Matrix I(int r, int c){
	assert(r==c);
	Matrix identity(r, c);
	for(int i = 0; i < r; i++){
		identity.matrix[i][i] = 1;
	}
	return identity;
}

double dotproduct(const Matrix A, const Matrix B){
	assert((A.rows==1) && (B.cols==1) && (A.cols==B.rows));
	double sum = 0;
	for(int i = 0; i < A.cols; i++){
		sum += A.matrix[0][i] * B.matrix[i][0];
	}
	return sum;
}

Matrix H(const Matrix u){
	// u is a column vector
	assert(u.cols==1);
	Matrix h = I(u.rows, u.rows) - (u * transpose(u)) * (2/dotproduct(transpose(u),u));
	return h;
}

// Householder Reflections, returns {Q, R} for square matrices
std::vector<Matrix> QR_decomposition(const Matrix A){
	assert(A.rows==A.cols);
	std::vector<Matrix> Qs;

	Matrix Aq = A;
	Matrix AI = I(A.cols, A.cols);
	double alpha;
	
	for(int i = 0; i < (A.cols-1); i++){
		Matrix Ap((A.rows-i),(A.cols-i));
		for(int j = 0; j < (A.cols-i); j++){
			for(int k = 0; k < (A.cols-i); k++){
				Ap.matrix[j][k] = Aq.matrix[j+i][k+i];
			}
		}
		// Form Householder matrix
		alpha = -sqrt(dotproduct(transpose(Ap.getcol(0)),Ap.getcol(0)));
		Matrix u = Ap.getcol(0);
		u.matrix[0][0] += alpha;
		Matrix Householder = H(u);
		
		Matrix Qi(A.cols, A.cols);
		// Pad Qi with identity matrix
		for(int l = 0; l < A.cols; l++){
			for(int m = 0; m < A.cols; m++){
				if((l < i) || (m < i)){
					Qi.matrix[l][m] = AI.matrix[l][m];
				}
				else{
					Qi.matrix[l][m] = Householder.matrix[l-i][m-i];
				} 
			}
		}
		
		Aq = Qi * Aq;
		
		Qs.push_back(Qi);
	}
	// Form Q, R
	Matrix Q = Qs[Qs.size()-1];
	for(int n = (Qs.size()-2); n >= 0; n--){
		Q = Qs[n]*Q;
	}
	
	// Form R
	Matrix R = A;
	for(int o = 0; o < Qs.size(); o++){
		R = Qs[o]*R;
	}

	std::vector<Matrix> result({Q, R});

	return result;
}

// Returns {D, Q} for square matrix
std::vector<Matrix> QR_diagonalize(const Matrix A, const double tol_add, const double tol_rat, const int maxiter){
	Matrix diag = A;
	Matrix temp;
	Matrix QT = I(diag.rows, diag.cols);
	std::vector<Matrix> QR(2);
	
	int count = 0;
	int good_count = 0;
	bool good = false;

	double nce_add = 2;
	double nce_rat = 2;
	while((!good) && (count < maxiter)){
		QR = QR_decomposition(diag);
		temp = transpose(QR[0]) * diag * QR[0];
		
		for(int i = 0; i < diag.rows; i++){
			nce_add = abs((temp.matrix[i][i]-diag.matrix[i][i])/((temp.matrix[i][i]+diag.matrix[i][i])/2));
			nce_rat = abs(temp.matrix[i][i]/diag.matrix[i][i]);
			std::cout << count << "   " << nce_add << '\n';
			std::cout << count << "   " << nce_rat << '\n';
			if((nce_add < tol_add) && (nce_rat < (1+tol_rat)) && (nce_rat > (1-tol_rat))){
				good_count++;
			}
		}
		std::cout << "good: " << good_count << '\n';
		if(good_count==diag.rows){
			good = true;
		}
		
		diag = temp;
		QT = transpose(QR[0]) * QT;
		good_count = 0;
		count+=1;
	}
	assert(count < maxiter);
	
	for(int j = 0; j < diag.rows; j++){
		for(int k = 0; k < diag.cols; k++){
			if(j!=k){
				diag.matrix[j][k] = 0;
			}
		}
	}

	Matrix Q = transpose(QT);

	std::vector<Matrix> result({diag, Q});
	return result;
}

Matrix inv_sqrt(const Matrix A){
	assert(A.rows==A.cols);
	std::vector<Matrix> QR = QR_diagonalize(A);
	for(int i = 0; i < QR[0].rows; i++){
		assert(QR[0].matrix[i][i] > 0);
		QR[0].matrix[i][i] = 1 / sqrt(QR[0].matrix[i][i]);
	}
	return QR[1] * QR[0] * transpose(QR[1]);
}
