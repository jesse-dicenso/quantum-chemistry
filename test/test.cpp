#include "../libint/2e.hpp"
#include "../libmol/mol.hpp"

using namespace std;

int main(){
	cout << setprecision(15);
	Molecule water("h2.inp");	
	
	Matrix Sij = overlap(water.AOs);
	Matrix Tij = kinetic(water.AOs);
	Matrix Vij = nuclear(water.AOs, water.Zvals, water.xyz);

	vector<vector<vector<vector<double>>>> eris = ERIs(water.AOs);
	
	Matrix Hcore = Tij + Vij;
	Matrix P = Hcore;
	Matrix g = G(P, eris);

	for(int i = 0; i < g.rows; i++){
		for(int j = 0; j < g.cols; j++){
			cout << "G" << i << j << " = " << g.matrix[i][j] << '\n';
		}
	}
	/*	
	for(int i = 0; i < eris.size(); i++){
		for(int j = 0; j < eris[i].size(); j++){
			for(int k = 0; k < eris[i][j].size(); k++){
				for(int l = 0; l < eris[i][j][k].size(); l++){
					cout << "g" << i << j << k << l << " = " << eris[i][j][k][l] << '\n';
				}
			}
		}
	}

	cout << "***********" << '\n';
	cout << "* OVERLAP *" << '\n';
	cout << "***********" << '\n';
	for(int i = 0; i < Sij.rows; i++){
		for(int j = 0; j < Sij.cols; j++){
			cout << "S" << i << j << " = " << Sij.matrix[i][j] << '\n';
		}
	}
	cout << '\n';
	cout << "***********" << '\n';
	cout << "* KINETIC *" << '\n';
	cout << "***********" << '\n';
	for(int i = 0; i < Tij.rows; i++){
		for(int j = 0; j < Tij.cols; j++){
			cout << "T" << i << j << " = " << Tij.matrix[i][j] << '\n';
		}
	}
	cout << '\n';
	cout << "***********" << '\n';
	cout << "* NUCLEAR *" << '\n';
	cout << "***********" << '\n';
	for(int i = 0; i < Tij.rows; i++){
		for(int j = 0; j < Tij.cols; j++){
			cout << "V" << i << j << " = " << Vij.matrix[i][j] << '\n';
		}
	}

	cout << '\n' << nucrepl(water.Zvals, water.xyz) << '\n';
	*/
}	
