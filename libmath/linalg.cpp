#include "linalg.hpp"

double& Matrix::operator()(int i, int j){
    if(!((i < rows) && (j < cols) && (i >= 0) && (j >= 0))){
        std::cout << "bad row: " << i << std::endl;
        std::cout << "bad col: " << j << std::endl;
        assert((i < rows) && (j < cols) && (i >= 0) && (j >= 0));
    }
    return data[cols * i + j];
}

const double& Matrix::operator()(int i, int j) const{
    if(!((i < rows) && (j < cols) && (i >= 0) && (j >= 0))){
        std::cout << "bad row: " << i << std::endl;
        std::cout << "bad col: " << j << std::endl;
        assert((i < rows) && (j < cols) && (i >= 0) && (j >= 0));
    }
    return data[cols * i + j];
}

Matrix Matrix::operator-() const{
	Matrix mat(rows, cols);
	for(int i = 0; i < rows * cols; i++){
		mat.data[i] = -data[i];
	}
	return mat;
}

Matrix& Matrix::operator=(const Matrix& A){
	rows = A.rows;
    cols = A.cols;
    data = A.data;
	return *this;
}

Matrix& Matrix::operator=(Matrix&& A){
	rows = A.rows;
    cols = A.cols;
    data = std::move(A.data);
	return *this;
}

Matrix Matrix::operator+(const Matrix& A) const{
	assert((rows==A.rows) && (cols==A.cols));
	Matrix sum(rows, cols);
	for(int i = 0; i < rows * cols; i++){
		sum.data[i] = data[i] + A.data[i];
	}
	return sum;
}

Matrix Matrix::operator-(const Matrix& A) const{
	assert((rows==A.rows) && (cols==A.cols));
	Matrix dif(rows, cols);
	for(int i = 0; i < rows * cols; i++){
		dif.data[i] = data[i] - A.data[i];
	}
	return dif;
}

Matrix Matrix::operator*(double c) const{
	Matrix mat(rows, cols);
	for(int i = 0; i < rows * cols; i++){
		mat.data[i] = data[i] * c;	
	}
	return mat;
}

Matrix Matrix::operator*(const Matrix& A) const{
	assert(cols==A.rows);
	Matrix product(rows, A.cols);
	for(int i = 0; i < product.rows; i++){
		for(int j = 0; j < product.cols; j++){
			for(int k = 0; k < cols; k++){
				product(i, j) += data[cols * i + k] * A(k, j);
			}
		}
	}
	return product;	
}

std::vector<double> Matrix::getRow(int r) const{
    assert((r < rows) && (r >= 0));
    std::vector<double> row(cols);
    for(int i = 0; i < cols; i++){
        row[i] = (*this)(r, i);
    }
    return row;
}

std::vector<double> Matrix::getCol(int c) const{
    assert((c < cols) && (c >= 0));
    std::vector<double> col(rows);
    for(int i = 0; i < rows; i++){
        col[i] = (*this)(i, c);
    }
    return col;
}

Matrix transpose(const Matrix& A){
	Matrix tp(A.cols, A.rows);
	for(int i = 0; i < A.rows; i++){
		for(int j = 0; j < A.cols; j++){
			tp(j, i) = A(i, j);
		}
	}
	return tp;
}

void mat_alloc(std::vector<Matrix>& AM, int size, int rows, int cols){
    for(int i = 0; i < size; i++){
        AM.emplace_back(rows, cols);
    }
}

Matrix I(int dim){
	Matrix identity(dim, dim);
	for(int i = 0; i < dim; i++){
		identity(i, i) = 1.0;
	}
	return identity;
}

Matrix zero(int r, int c){
	Matrix z(r, c);
	return z;
}

double Tr(const Matrix& A){
	assert(A.rows==A.cols);
	double sum = 0;
	for(int i = 0; i < A.rows; i++){
		sum += A(i, i);
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
			a[i+n*j] = A(i, j);
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
	for(int i = 0; i < n; i++){
		diag(i, i) = w[i];
	}
	
	Matrix Q(n, n);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			Q(i, j) = a[i+n*j];
		}
	}	

	std::vector<Matrix> result({diag, Q});
	return result;
}

Matrix m_sqrt(const Matrix& A){
	assert(A.rows==A.cols);
	std::vector<Matrix> QR = diagonalize(A);
	for(int i = 0; i < QR[0].rows; i++){
		assert(QR[0](i, i) > 0);
		QR[0](i, i) = sqrt(QR[0](i, i));
	}
	return QR[1] * QR[0] * transpose(QR[1]);
}

Matrix m_inv_sqrt(const Matrix& A){
	assert(A.rows==A.cols);
	std::vector<Matrix> QR = diagonalize(A);
	for(int i = 0; i < QR[0].rows; i++){
		assert(QR[0](i, i) > 0);
		QR[0](i, i) = 1 / sqrt(QR[0](i, i));
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
			a[i+n*j] = A(i, j);
		}
	}
	int lda = n;				// leading dimension of A
	std::vector<int> ipiv(n);	// Details of D
	std::vector<double> b(n);	// B matrix, column-major
	for(int i = 0; i < n; i++){
		b[i] = B(i, 0);
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

/*void Matrix::printMatrix(int w){
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			std::cout << std::setw(w) << data(i * rows + j) << std::setw(w);
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}*/
