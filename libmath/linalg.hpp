#ifndef LINALGHEADERDEF
#define LINALGHEADERDEF

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

// LAPACK function for computing eigenvalues/eigenvectors
extern "C" {
	void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO);
	void dsysv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* WORK, int* LWORK, int* INFO);
}

class Matrix{
	public:
		Matrix(int r, int c, bool sym=false);
		Matrix(const Matrix& A);
		Matrix();
		~Matrix();

		int rows;
		int cols;
		double **matrix;
		bool isSymmetric;
		
		void printMatrix(int w=20);
		Matrix getrow(int i);
		Matrix getcol(int i);

		Matrix operator-() const;
		Matrix& operator=(const Matrix&A);
		Matrix operator+(const Matrix&A) const;
		Matrix operator-(const Matrix& A) const;
		Matrix operator*(double c) const;
		Matrix operator*(const Matrix& A) const;
};

Matrix I(int r, int c);
Matrix zero(int r, int c);

Matrix transpose(const Matrix A);

double Tr(Matrix A);

std::vector<Matrix> diagonalize(Matrix A);
Matrix m_sqrt(const Matrix A);
Matrix m_inv_sqrt(const Matrix A);

std::vector<double> sym_linear_solve(Matrix A, Matrix B, int* icd);

#endif
