#ifndef MOLHEADERDEF
#define MOLHEADERDEF

#include "../libgf/gf.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

class Molecule{
	public:
		Molecule(const char *file);
		~Molecule();

		int Natoms;
		int charge;
		int* Zvals;
		int Nelec;
		double** xyz;
		GF** MO;
};

#endif
