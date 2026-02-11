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
		Matrix(int rows_, int cols_) : rows(rows_), cols(cols_), data(rows_ * cols_) {}
		Matrix(const Matrix& A) : rows(A.rows), cols(A.cols), data(A.data) {}
		Matrix() : rows(0), cols(0) {}

		int rows;
		int cols;
        std::vector<double> data;
		
        double& operator()(int i, int j);
        const double& operator()(int i, int j) const;
		Matrix operator-() const;
		Matrix& operator=(const Matrix& A);
		Matrix& operator=(Matrix&& A);
		Matrix operator+(const Matrix& A) const;
		Matrix operator-(const Matrix& A) const;
		Matrix operator*(double c) const;
		Matrix operator*(const Matrix& A) const;
        
        std::vector<double> getRow(int r) const;
        std::vector<double> getCol(int c) const;
        /*Matrix getRow(int r) const;
        Matrix getCol(int c) const;*/

        void resize(int r_new, int c_new);
};

Matrix I(int r, int c);
Matrix zero(int r, int c);
void mat_alloc(std::vector<Matrix>& AM, int size, int rows, int cols);
Matrix transpose(const Matrix& A);
double Tr(const Matrix& A);
std::vector<Matrix> diagonalize(const Matrix& A);
Matrix m_sqrt(const Matrix& A);
Matrix m_inv_sqrt(const Matrix& A);
std::vector<double> sym_linear_solve(const Matrix& A, const Matrix& B, int* icd);

class Tensor3{
    public:
        Tensor3(int dim1_, int dim2_, int dim3_) : dim1(dim1_), dim2(dim2_), dim3(dim3_), data(dim1_ * dim2_ * dim3_) {}
        Tensor3(const Tensor3& A): dim1(A.dim1), dim2(A.dim2), dim3(A.dim3), data(A.data) {}
        Tensor3() : dim1(0), dim2(0), dim3(0) {}

        int dim1;
        int dim2;
        int dim3;
        std::vector<double> data;

        double& operator()(int i, int j, int k);
        const double& operator()(int i, int j, int k) const;
        Tensor3 operator-() const;
        Tensor3& operator=(const Tensor3& A);
        Tensor3& operator=(Tensor3&& A);
        Tensor3 operator+(const Tensor3& A) const;
        Tensor3 operator-(const Tensor3& A) const;
        Tensor3 operator*(double c) const;
        
        void resize(int dim1_new, int dim2_new, int dim3_new);
};

#endif
