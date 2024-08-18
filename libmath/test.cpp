#include "linalg.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

int main(){
	cout << fixed;
	cout << setprecision(15);
	Matrix F(2,2);
	F.matrix[0][0] = -0.3655;
	F.matrix[1][0] = -0.5939;
	F.matrix[0][1] = -0.5939;
	F.matrix[1][1] = -0.3655;
	
	double st = 1.0 / sqrt(2*(1+0.6593));
	double sm = 1.0 / sqrt(2*(1-0.6593));

	Matrix C(2,2);
	C.matrix[0][0] = st;
	C.matrix[1][0] = st;
	C.matrix[0][1] = sm;
	C.matrix[1][1] = -sm;
	C.printMatrix();

	Matrix S(2,2);
	S.matrix[0][0] = 1;	
	S.matrix[1][0] = 0.6593;
	S.matrix[0][1] = 0.6593;
	S.matrix[1][1] = 1;

	Matrix e(2,2);
	C.matrix[0][0] = -1.2528;
	C.matrix[1][0] = 0;
	C.matrix[0][1] = 0;
	C.matrix[1][1] = -0.4756;

	(F*C-S*C*e).printMatrix();

}
