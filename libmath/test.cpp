#include "linalg.hpp"
#include <iostream>

using namespace std;

int main(){
	Matrix A(3,3);
	A.matrix[0][0] = 12;
	A.matrix[1][0] = 6;
	A.matrix[2][0] = -4;
	A.matrix[0][1] = -51;
	A.matrix[1][1] = 167;
	A.matrix[2][1] = 24;
	A.matrix[0][2] = 4;
	A.matrix[1][2] = -68;
	A.matrix[2][2] = -41;

	vector<Matrix> QR = QR_decomposition(A);
	cout << "Q : \n";
	QR[0].printMatrix();
	cout << '\n';
	cout << "R : \n";
	QR[1].printMatrix();
	cout << '\n';
}
