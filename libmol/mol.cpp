#include "mol.hpp"

Molecule::Molecule(const char *file){
	std::ifstream read(file);
	assert(read.good());
	
	read >> Natoms >> charge;
	Nelec = 0 - charge;
	Zvals = new int[Natoms];
	xyz  = new double*[Natoms];
	for(int i = 0; i < Natoms; i++){
		xyz[i] = new double[3];
	}
	for(int j = 0; j < Natoms; j++){
		read >> Zvals[j] >> xyz[j][0] >> xyz[j][1] >> xyz[j][2];
		Nelec += Zvals[j];
	}
	// closed-shell
	assert((Nelec%2)==0);

	std::string basis;
	read >> basis;
	// only STO-3G available...sorry
	assert(basis=="STO-3G");
	read.close();
}

Molecule::~Molecule(){
	for(int i = 0; i < Natoms; i++){
		delete[] xyz[i];
	}
	delete[] xyz;
	delete[] Zvals;
	
}
