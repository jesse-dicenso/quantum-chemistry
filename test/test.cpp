#include "../libint/1e.hpp"

using namespace std;

int main(){
	vector<double> a({0.168856,0.623913,3.42525});
	vector<double> d({0.444635,0.535328,0.154329});
	vector<double> ctr1({0.0, 1.4194772274, -0.87926122007});
	vector<double> ctr2({0.0, -1.4194772274, -0.87926122007});
	//vector<double> ctr1({0.0, 1.4, 0.0});
	//vector<double> ctr2({0.0, 0.0, 0.0});
	
	vector<int> L({0,0,0});
	PGF p1(a[0],ctr1,L);
	PGF p2(a[1],ctr1,L);
	PGF p3(a[2],ctr1,L);
	vector<PGF> p({p1, p2, p3});
	PGF q1(a[0],ctr2,L);
	PGF q2(a[1],ctr2,L);
	PGF q3(a[2],ctr2,L);
	vector<PGF> q({q1, q2, q3});
	CGF c1(p, d);
	CGF c2(q, d);
	
	//cout << "Spg = " << setprecision(15) << S(p1, p1) << '\n'; 
	cout << "S11 = " << setprecision(15) << S(c1, c1) << '\n'; 
	cout << "S12 = " << setprecision(15) << S(c1, c2) << '\n'; 
	cout << "S21 = " << setprecision(15) << S(c2, c1) << '\n'; 
	cout << "S22 = " << setprecision(15) << S(c2, c2) << '\n'; 
}	
