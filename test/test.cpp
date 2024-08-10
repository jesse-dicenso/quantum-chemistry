#include "../libint/1e.hpp"

using namespace std;

int main(){
	vector<double> a({0.168856,0.623913,3.42525});
	vector<double> d({0.444635,0.535328,0.154329});
	// Water
	vector<double> ctr1({0.0, 1.4194772274, -0.87926122007});	// H
	vector<double> ctr3({0.0, 0.0, 0.21981483259});			// O
	vector<double> ctr2({0.0, -1.4194772274, -0.87926122007});	// H
	// H2
	//vector<double> ctr1({0.0, 1.4, 0.0});
	//vector<double> ctr2({0.0, 0.0, 0.0});	
	vector<int> Ls({0,0,0});

	vector<double> v1({1.0,2.0,3.0});
	vector<double> v2({1.0,1.0,1.0});

	GF g1(a,d,ctr1,Ls);
	GF g2(a,d,ctr2,Ls);
	
	cout << "S11 = " << setprecision(15) << S(g1, g1) << '\n'; 
	cout << "S12 = " << setprecision(15) << S(g1, g2) << '\n'; 
	cout << "S21 = " << setprecision(15) << S(g2, g1) << '\n'; 
	cout << "S22 = " << setprecision(15) << S(g2, g2) << '\n';
	cout << "T11 = " << setprecision(15) << T(g1, g1) << '\n'; 
	cout << "T12 = " << setprecision(15) << T(g1, g2) << '\n'; 
	cout << "T21 = " << setprecision(15) << T(g2, g1) << '\n'; 
	cout << "T22 = " << setprecision(15) << T(g2, g2) << '\n';
}	
