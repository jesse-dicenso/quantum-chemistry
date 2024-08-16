#include "linalg.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

int main(){
	cout << fixed;
	cout << setprecision(15);
	Matrix A(2,2);
	A.matrix[0][0] = -2.459500456944964;
	A.matrix[1][0] = 0.011616567307498;
	A.matrix[0][1] = 0.011616567307498;
	A.matrix[1][1] = -2.459506040336783;
	/*
	Matrix A(2,2);
	A.matrix[0][0] = 4;
	A.matrix[1][0] = 1;
	A.matrix[0][1] = -1;
	A.matrix[1][1] = 3;
	*/
	vector<Matrix> D = QR_diagonalize(A);
	
	cout << "D : \n";
	D[0].printMatrix();
	cout << '\n';
	cout << "Q : \n";
	D[1].printMatrix();
	cout << '\n';
}
