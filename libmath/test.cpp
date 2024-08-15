#include "linalg.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

int main(){
	cout << setprecision(15);
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

	vector<Matrix> QR = QR_diagonalize(A);
	cout << "D : \n";
	QR[0].printMatrix();
	cout << '\n';
	cout << "Q : \n";
	QR[1].printMatrix();
	cout << '\n';
}
