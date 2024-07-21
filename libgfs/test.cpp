#include <cmath>
#include <iostream>
#include <iomanip>
#include "gf.hpp"
#include <vector>

using namespace std;

int main(){
	vector<double> center({1.0994558, 2.598574, -12.000000123});
	PGF p1(1.27, center[0], center[1], center[2], 1, 0, 0);
	PGF p2(1.0, center[0], center[1], center[2], 1, 0, 0);
	PGF p3(1.35, center[0], center[1], center[2], 1, 0, 0);
	vector<PGF> pv({p1, p2, p3});
	vector<double> dv({0.67978, 0.012334, 0.02221});
	CGF c1(3, pv, dv);
	c1.printcgf();
}
