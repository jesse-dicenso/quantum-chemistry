#ifndef MOLHEADERDEF
#define MOLHEADERDEF

#include "../libgf/gf.hpp"

#include <fstream>
#include <iostream>
#include <string>

class Molecule{
	public:
		Molecule(const char *file);

		int Natoms;
		std::vector<int> Zvals;
		int Nelec;
		std::vector<std::vector<double>> xyz;
		std::vector<GF> AOs;
};

std::vector<GF> AOfunctions(std::string basis, int Zval, std::vector<double> xyz);

#endif
