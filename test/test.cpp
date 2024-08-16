#include "../libint/2e.hpp"
#include "../libmol/mol.hpp"

using namespace std;

int main(){
	cout << setprecision(3);
	cout << fixed;
	int pr = 10;
	Molecule water("H2O.inp");	
	
	Matrix Sij = overlap(water.AOs);
	Matrix X = inv_sqrt(Sij);

	Matrix Tij = kinetic(water.AOs);
	Matrix Vij = nuclear(water.AOs, water.Zvals, water.xyz);

	vector<vector<vector<vector<double>>>> eris = ERIs(water.AOs);

	cout << "S\n";
	Sij.printMatrix(pr);

	cout << "T\n";
	Tij.printMatrix(pr);

	cout << "V\n";
	Vij.printMatrix(pr);
		
	Matrix Hcore = Tij + Vij;

	cout << "\nHcore\n";
	Hcore.printMatrix(pr);
	//Matrix g = G(P, eris);
}	
