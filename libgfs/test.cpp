#include <cmath>
#include <iostream>
#include <iomanip>
#include "gf.hpp"
//#include "../libint/0e.hpp"
#include <vector>

using namespace std;

int main(){
	vector<double> center1({1.5994558, 2.598574, 12.300000123});
	vector<double> center2({1.84558, 2.498574, 12.000000123});
	PGF p1(1.27, center1, 2, 3, 4);
	PGF p2(1.34, center1, 2, 3, 4);
	PGF p3(1.05, center1, 2, 3, 4);
	vector<PGF> pg1({p1, p2, p3});
	PGF q1(1.52, center2, 1, 2, 1);
	PGF q2(1.67, center2, 1, 2, 1);
	PGF q3(1.43, center2, 1, 2, 1);
	vector<PGF> pg2({q1, q2, q3});
	vector<double> d1({0.32333, 0.29349, 0.9111});
	vector<double> d2({0.7634, 0.28855, 0.1993443});
	vector<double> d3({0.32333, 0.293491, 0.9111});
	CGF c1(3, pg1, d1);
	CGF c2(3, pg2, d2);
	CGF c3(3, pg1, d3);
	cout << c1.getN() <<'\n';
	cout << c2.getN() << '\n';
	/*
	cout << setprecision(15) << overlap(p1, p1) << '\n';
	cout << setprecision(15) << overlap(p1, p2) << '\n';
	cout << setprecision(15) << overlap(p1, p3) << '\n';
	cout << setprecision(15) << overlap(p1, q1) << '\n';
	cout << setprecision(15) << overlap(p1, q2) << '\n';
	cout << setprecision(15) << overlap(p1, q3) << '\n';
	*/
	cout << setprecision(15) << overlap(c1, c1) << '\n';
	cout << setprecision(15) << overlap(c1, c2) << '\n';
	cout << setprecision(15) << overlap(c2, c1) << '\n';
	cout << setprecision(15) << overlap(c2, c2) << '\n';
	cout << setprecision(15) << overlap(c1, c3) << '\n';
	/*
	vector<double> center({0.0, 0.0, 0.0});
	PGF t1(1.0, center, 2, 1, 3);
	PGF t2(1.0, center, 2, 1, 3);
	vector<PGF> ts({t1, t2});
	vector<double> ds({0.5, 0.5});
	CGF t12(2, ts, ds);
	cout << t1.getN() << '\n';
	cout << t12.getN() << '\n';
	*/
}	
