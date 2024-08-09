#include "../libint/1e.hpp"

using namespace std;

int main(){
	vector<double> a({0.168856,0.623913,3.42525});
	vector<double> d({0.444635,0.535328,0.154329});
	vector<double> ctr1({0, 1.4194772274, -0.87926122007});
	vector<double> ctr2({0, -1.4194772274, -0.87926122007});
	//vector<double> ctr1({0.0, 1.4, 0.0});
	//vector<double> ctr2({0.0, 0.0, 0.0});	
	vector<int> L({0,0,0});

	GF g1(a,d,ctr1,L);
	GF g2(a,d,ctr2,L);
	
	cout << "S11 = " << setprecision(15) << S(g1, g1) << '\n'; 
	cout << "S12 = " << setprecision(15) << S(g1, g2) << '\n'; 
	cout << "S21 = " << setprecision(15) << S(g2, g1) << '\n'; 
	cout << "S22 = " << setprecision(15) << S(g2, g2) << '\n'; 
}	
