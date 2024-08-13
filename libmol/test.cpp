#include "mol.hpp"

using namespace std;

int main(){
	Molecule m("input");
	cout << m.Natoms << " " << m.charge << '\n';
	for(int i = 0; i < m.Natoms; i++){
		cout << m.Zvals[i] << "  " << m.xyz[i][0] << "   " << m.xyz[i][1] << "   " << m.xyz[i][2] << '\n';
	}
}
