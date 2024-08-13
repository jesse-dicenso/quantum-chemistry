#include "mol.hpp"

Molecule::Molecule(const char *file){
	std::ifstream read(file);
	assert(read.good());
	
	read >> Natoms >> charge;
	Zvals = new int[Natoms];
	xyz  = new double*[Natoms];
	for(int i = 0; i < Natoms; i++){
		xyz[i] = new double[3];
	}
	for(int j = 0; j < Natoms; j++){
		read >> Zvals[j] >> xyz[j][0] >> xyz[j][1] >> xyz[j][2];
	}
	read.close();
}

Molecule::~Molecule(){
	for(int i = 0; i < Natoms; i++){
		delete[] xyz[i];
	}
	delete[] xyz;
	delete[] Zvals;
	
}
