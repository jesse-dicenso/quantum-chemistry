#include "main.hpp"

using namespace std;

int main(){
	// GENERAL ELECTRONIC STRUCTURE PROGRAM
	cout << setprecision(15);
	Molecule M("INPUT");
	int K = M.AOs.size();
	
	// Nuclear repulsion energy	
	double nuc = nucrepl(M.Zvals, M.xyz);

	// One electron integrals
	Matrix s = overlap(M.AOs);
	Matrix t = kinetic(M.AOs);
	Matrix v = nuclear(M.AOs, M.Zvals, M.xyz);

	Matrix hcore = t + v;
	
	// Two electron integrals
	vector<vector<vector<vector<double>>>> g = ERIs(M.AOs);	
}
