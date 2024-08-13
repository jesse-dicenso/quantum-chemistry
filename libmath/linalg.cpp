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

Matrix::~Matrix(){
	for(int i = 0; i < rows; i++){
		delete[] matrix[i];
	}
	delete[] matrix;
}

void Matrix::printMatrix(){
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			std::cout << matrix[i][j] << "    ";
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
	assert((rows==A.rows) && (cols==A.cols));
	isSymmetric = A.isSymmetric;
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			matrix[i][j] == A.matrix[i][j];
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
