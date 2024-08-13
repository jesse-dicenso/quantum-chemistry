#include "../libint/2e.hpp"

using namespace std;

int main(){
	vector<double> a({0.168856,0.623913,3.42525});
	vector<double> d({0.444635,0.535328,0.154329});
	// Water
	//vector<double> ctr1({0.0, 1.4194772274, -0.87926122007});	// H
	//vector<double> ctr3({0.0, 0.0, 0.21981483259});			// O
	//vector<double> ctr2({0.0, -1.4194772274, -0.87926122007});	// H
	// H2
	vector<double> ctr1({0.0, 1.4, 0.0});
	vector<double> ctr2({0.0, 0.0, 0.0});	
	vector<int> Ls({0,0,0});

	GF g1(a,d,ctr1,Ls);
	GF g2(a,d,ctr2,Ls);
	
	cout << "S = [ " << S(g1, g1); 
	cout <<      " " << S(g1, g2) << " ]\n"; 
	cout << "    [ " << S(g2, g1); 
	cout <<      " " << S(g2, g2) << " ]\n\n";
	
	cout << "T = [ " << T(g1, g1); 
	cout <<      " " << T(g1, g2) << " ]\n"; 
	cout << "    [ " << T(g2, g1); 
	cout <<      " " << T(g2, g2) << " ]\n\n";
	
	cout << "V111 = -" << V(g1, g1, ctr1) << '\n';
	cout << "V211 = -" << V(g2, g1, ctr1) << '\n';
	cout << "V121 = -" << V(g1, g2, ctr1) << '\n';
	cout << "V221 = -" << V(g2, g2, ctr1) << '\n';
	cout << "V112 = -" << V(g1, g1, ctr2) << '\n';
	cout << "V112 = -" << V(g2, g1, ctr2) << '\n';
	cout << "V112 = -" << V(g1, g2, ctr2) << '\n';
	cout << "V112 = -" << V(g2, g2, ctr2) << "\n\n";
	
	cout << "g1111 = " << G(g1, g1, g1, g1) << '\n';
	cout << "g1122 = " << G(g1, g1, g2, g2) << '\n';
	cout << "g2222 = " << G(g2, g2, g2, g2) << '\n';
}	
