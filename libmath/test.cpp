#include "linalg.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

int main(){
	Matrix A(2,2);
	A.matrix[0][0] = -2.45950324182410;
	A.matrix[1][0] = 1.0116165673111346;
	A.matrix[0][1] = -1.0116165673111346;
	A.matrix[1][1] = -2.45950325545770;

	vector<Matrix> D = QR_diagonalize(A, 0.001, 0.001);
	
	cout << "D : \n";
	D[0].printMatrix();
	cout << '\n';
	cout << "Q : \n";
	D[1].printMatrix();
	cout << '\n';
}
