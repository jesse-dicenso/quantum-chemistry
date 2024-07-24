#include <cmath>
#include <iostream>
#include <iomanip>
#include "gf.hpp"
#include <vector>

using namespace std;

int main(){
	/*
	vector<double> center({1.0994558, 2.598574, -12.000000123});
	PGF p1(1.27, center[0], center[1], center[2], 2, 1, 0);
	PGF p2(1.0, center[0], center[1], center[2], 2, 1, 0);
	PGF p3(1.35, center[0], center[1], center[2], 2, 1, 0);
	vector<PGF> pv({p1, p2, p3});
	vector<double> dv({0.67978, 0.012334, 0.02221});
	p1.printpgf();
	CGF c1(3, pv, dv);
	c1.printcgf();
	*/
	vector<double> xyz1({1.5, 2.5, 3.5});
	vector<double> xyz2({-0.5, 0.5, -3.0});
	vector<double> P(3);
	PGF p4(1.38, xyz1, 2, 1, 0);
	PGF p5(1.27, xyz2, 1, 1, 1);
	PGF p6(1.38, xyz1, 2, 1, 0);
	PGF p7(1.38, xyz2, 2, 1, 0);
	PGF p8(1.38, xyz1, 2, 2, 0);

	cout << (p4==p5) << '\n';
	cout << (p4==p6) << '\n';
	cout << (p4==p7) << '\n';
	cout << (p4==p8) << '\n';
}	
