#include "linalg.hpp"
#include <iostream>

using namespace std;

int main(){
	Matrix A(3,3);
	A.matrix[0][0] = 1;
	A.matrix[1][0] = -2;
	A.matrix[2][0] = 3;
	A.matrix[0][1] = -2;
	A.matrix[1][1] = 1;
	A.matrix[2][1] = 4.5;
	A.matrix[0][2] = 3;
	A.matrix[1][2] = 4.5;
	A.matrix[2][2] = 1;

	Matrix B(3,1);
	B.matrix[0][0] = 1;
	B.matrix[1][0] = 1;
	B.matrix[2][0] = 1;

	vector<double> sol = sym_linear_solve(A, B);

	for(int i = 0; i < 3; i++){
		cout << sol[i] << '\n';
	}
}
