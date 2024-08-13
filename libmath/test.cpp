#include "linalg.hpp"
#include <iostream>

using namespace std;

int main(){
	Matrix A(3,3);
	A.matrix[0][0] = 2.5;
	A.matrix[1][1] = -12;
	A.matrix[2][2] = 0.5;
	A.matrix[2][1] = 1;
	Matrix B(3,3);
	B.matrix[0][1] = 1;
	B.matrix[1][1] = 1;
	B.matrix[2][2] = 1;
	(A*B).printMatrix();
}
