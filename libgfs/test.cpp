#include <cmath>
#include <iostream>
#include <iomanip>
#include "gf.hpp"
#include <vector>

using namespace std;

int main(){
	vector<double> center1({1.0994558, 2.598574, -12.000000123});
	vector<double> center2({-1.0994558, -2.598574, 12.000000123});
	PGF p1(1.27, center1, 2, 1, 0);
	PGF p2(1.27, center2, 1, 1, 1);
	cout << "K = " << K(p1, p2) << '\n';
	cout << "P = " << P(p1, p2)[0] << "   " << P(p1, p2)[1] << "   " << P(p1, p2)[2] << '\n';
}	
