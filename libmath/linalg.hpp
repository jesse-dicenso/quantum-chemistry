#ifndef LINALGHEADERDEF
#define LINALGHEADERDEF

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

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
double dot(const Matrix A, const Matrix B);

// QR algorithm via Householder Reflections
Matrix H(const Matrix u);
std::vector<Matrix> QR_decomposition(const Matrix A);
std::vector<Matrix> QR_diagonalize(const Matrix A, const double tol=1e-8, const int maxiter=500);

Matrix inv_sqrt(const Matrix A);

#endif
