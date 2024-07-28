#include "../libint/1e.hpp"

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
	cout << "S(p1, p1) = " << S(p1, p1) << '\n';
	cout << "S(p1, p2) = " << S(p1, p2) << '\n';
	cout << "S(p2, p1) = " << S(p2, p1) << '\n';
	cout << "S(p1, q1) = " << S(p1, q1) << '\n';
	cout << "S(q1, p1) = " << S(q1, p1) << '\n';
	cout << "\nS(c1, c1) = " << S(c1, c1) << '\n';
	cout << "S(c1, c2) = " << S(c1, c2) << '\n';	
	cout << "S(c2, c1) = " << S(c2, c1) << '\n';

	cout << "\nT(p1, p1) = " << T(p1, p1) << '\n';
	cout << "T(p1, p2) = " << T(p1, p2) << '\n';
	cout << "T(p2, p1) = " << T(p2, p1) << '\n';
	cout << "T(p1, q1) = " << T(p1, q1) << '\n';
	cout << "T(q1, p1) = " << T(q1, p1) << '\n';
	cout << "\nT(c1, c1) = " << T(c1, c1) << '\n';
	cout << "T(c1, c2) = " << T(c1, c2) << '\n';	
	cout << "T(c2, c1) = " << T(c2, c1) << '\n';
}	
