#ifndef MOLHEADERDEF
#define MOLHEADERDEF

#include <cassert>
#include <fstream>
#include <iostream>

class Molecule{
	public:
		Molecule(const char *file);
		~Molecule();

		int Natoms;
		int charge;
		int* Zvals;
		double** xyz;
};

#endif
