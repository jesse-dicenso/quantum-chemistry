#include "../libint/buildmatrix.hpp"

using namespace std;

int main(){
	cout << setprecision(15);
	// STO-3G //
	// * H *
	// 1s
	vector<double> aH1s({0.168856,0.623913,3.42525});
	vector<double> dH1s({0.444635,0.535328,0.154329});
	// * O *
	// 1s
	vector<double> aO1s({0.1307093214e3,0.2380886605e2,0.6443608313e1});
	vector<double> dO1s({0.1543289673,0.5353281423,0.4446345422});
	// 2s
	vector<double> aO2sp({0.5033151319e1,0.1169596125e1,0.3803889600});
	vector<double> dO2s ({-0.9996722919e-1,0.3995128261,0.7001154689});
	// 2p	
	vector<double> dO2p({0.1559162750,0.6076837186,0.3919573931});
	// Water
	vector<double> ctr1({0.0, 1.4194772274, -0.87926122007});	// H
	vector<double> ctr2({0.0, -1.4194772274, -0.87926122007});	// H
	vector<double> ctr3({0.0, 0.0, 0.21981483259});			// O
	// H2
	//vector<double> ctr1({0.0, 1.4, 0.0});
	//vector<double> ctr2({0.0, 0.0, 0.0});	
	vector<int> Ls ({0,0,0});
	vector<int> Lpx({1,0,0});
	vector<int> Lpy({0,1,0});
	vector<int> Lpz({0,0,1});
	
	vector<int> Z({8,1,1});
	vector<vector<double>> pos({ctr3,ctr1,ctr2});

	GF g1(aO1s,dO1s,ctr3,Ls);
	GF g2(aO2sp,dO2s,ctr3,Ls);
	GF g3(aO2sp,dO2p,ctr3,Lpx);
	GF g4(aO2sp,dO2p,ctr3,Lpy);
	GF g5(aO2sp,dO2p,ctr3,Lpz);
	GF g6(aH1s,dH1s,ctr1,Ls);
	GF g7(aH1s,dH1s,ctr2,Ls);

	vector<GF> gfs({g1,g2,g3,g4,g5,g6,g7});
	
	Matrix Sij = overlap(gfs);
	Matrix Tij = kinetic(gfs);
	Matrix Vij = nuclear(gfs, Z, pos);

	cout << "***********" << '\n';
	cout << "* OVERLAP *" << '\n';
	cout << "***********" << '\n';
	for(int i = 0; i < Sij.rows; i++){
		for(int j = 0; j < Sij.cols; j++){
			cout << "S" << i << j << " = " << Sij.matrix[i][j] << '\n';
		}
	}
	cout << "***********" << '\n';
	cout << "* KINETIC *" << '\n';
	cout << "***********" << '\n';
	for(int i = 0; i < Tij.rows; i++){
		for(int j = 0; j < Tij.cols; j++){
			cout << "T" << i << j << " = " << Tij.matrix[i][j] << '\n';
		}
	}
	cout << "***********" << '\n';
	cout << "* NUCLEAR *" << '\n';
	cout << "***********" << '\n';
	for(int i = 0; i < Tij.rows; i++){
		for(int j = 0; j < Tij.cols; j++){
			cout << "V" << i << j << " = " << Vij.matrix[i][j] << '\n';
		}
	}
}	
