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

void Matrix::printMatrix(int w){
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			std::cout << std::setw(w) << matrix[i][j] << std::setw(w);
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

Matrix Matrix::operator-() const{
	Matrix mat(rows, cols);
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			mat.matrix[i][j] = -matrix[i][j];
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

Matrix transpose(const Matrix& A){
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
		identity.matrix[i][i] = 1.0;
	}
	return identity;
}

Matrix zero(int r, int c){
	Matrix z(r, c);
	for(int i = 0; i < r; i++){
		for(int j = 0; j < c; j++){
			z.matrix[i][i] = 0.0;
		}
	}
	return z;
}

double Tr(const Matrix& A){
	assert(A.rows==A.cols);
	double sum = 0;
	for(int i = 0; i < A.rows; i++){
		sum += A.matrix[i][i];
	}
	return sum;
}

std::vector<Matrix> diagonalize(const Matrix& A){
	assert(A.rows==A.cols);
	// * DSYEV options * //
	char jobz = 'V';		// eigenvalues & eigenvectors
	char uplo = 'U';		// upper triangular
	int n = A.cols;			// matrix size	
	std::vector<double> a(n*n);	// matrix, column-major
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			a[i+n*j] = A.matrix[i][j];
		}
	}
	int lda = n;			// leading dimension
	std::vector<double> w(n);	// eigenvalues
	int info;			// successful if info == 0
	
	// workspace query
	double works;
	int lwork = -1;
	dsyev_(&jobz, &uplo, &n, a.data(), &lda, w.data(), &works, &lwork, &info);
	
	if(info != 0){
		std::cerr << "Workspace query err, info = " << info << '\n';
		return {};
	}

	lwork = (int)works;
	if(lwork <= 0){
		std::cerr << "Illegal lwork value, info = " << info << '\n';
		return {};
	}

	assert(lwork>0);
	std::vector<double> work(lwork);

	// Calculate eigenvalues/vectors
	dsyev_(&jobz, &uplo, &n, a.data(), &lda, w.data(), work.data(), &lwork, &info);	
	if(info != 0){
		std::cerr << "info != 0\n";
		return {};
	}

	// Return output as Matrix types
	Matrix diag(n, n);
	for(int i = 0; i < diag.rows; i++){
		diag.matrix[i][i] = w[i];
	}
	
	Matrix Q(n, n);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			Q.matrix[i][j] = a[i+n*j];
		}
	}	

	std::vector<Matrix> result({diag, Q});
	return result;
}

Matrix m_sqrt(const Matrix& A){
	assert(A.rows==A.cols);
	std::vector<Matrix> QR = diagonalize(A);
	for(int i = 0; i < QR[0].rows; i++){
		assert(QR[0].matrix[i][i] > 0);
		QR[0].matrix[i][i] = sqrt(QR[0].matrix[i][i]);
	}
	return QR[1] * QR[0] * transpose(QR[1]);
}

Matrix m_inv_sqrt(const Matrix& A){
	assert(A.rows==A.cols);
	std::vector<Matrix> QR = diagonalize(A);
	for(int i = 0; i < QR[0].rows; i++){
		assert(QR[0].matrix[i][i] > 0);
		QR[0].matrix[i][i] = 1 / sqrt(QR[0].matrix[i][i]);
	}
	return QR[1] * QR[0] * transpose(QR[1]);
}

std::vector<double> sym_linear_solve(const Matrix& A, const Matrix& B, int* icd){
	assert(A.rows==A.cols);
	assert((B.rows==A.rows) && (B.cols==1));
	// * DSYSV options * //
	char uplo = 'U';			// upper triangular
	int n = A.cols;				// A matrix size
	int nrhs = 1;				// 1 column B matrix
	std::vector<double> a(n*n);	// A matrix, column-major
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			a[i+n*j] = A.matrix[i][j];
		}
	}
	int lda = n;				// leading dimension of A
	std::vector<int> ipiv(n);	// Details of D
	std::vector<double> b(n);	// B matrix, column-major
	for(int i = 0; i < n; i++){
		b[i] = B.matrix[i][0];
	}
	int ldb = n;				// leading dimension of B
	int info;					// successful if info == 0
	
	// workspace query
	double works;
	int lwork = -1;
	dsysv_(&uplo, &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, &works, &lwork, &info);
	
	if(info != 0){
		std::cerr << "Workspace query err, info = " << info << '\n';
		return {};
	}

	lwork = (int)works;
	if(lwork <= 0){
		std::cerr << "Illegal lwork value, info = " << info << '\n';
		return {};
	}

	assert(lwork>0);
	std::vector<double> work(lwork);

	// Solve AX = B
	dsysv_(&uplo, &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, work.data(), &lwork, &info);
	*icd = info;

	return b;
}
