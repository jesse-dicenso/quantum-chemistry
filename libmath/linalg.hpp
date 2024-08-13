#ifndef LINALGHEADERDEF
#define LINALGHEADERDEF

#include <cassert>
#include <cmath>
#include <iostream>

class Matrix{
	public:
		Matrix(int r, int c, bool sym=false);
		Matrix(const Matrix& A);
		~Matrix();

		int rows;
		int cols;
		double **matrix;
		bool isSymmetric;
		
		void printMatrix();		

		Matrix operator-() const;
		Matrix& operator=(const Matrix&A);
		Matrix operator+(const Matrix&A) const;
		Matrix operator-(const Matrix& A) const;
		Matrix operator*(double c) const;
		Matrix operator*(const Matrix& A) const;
};

Matrix transpose(const Matrix A);

#endif
