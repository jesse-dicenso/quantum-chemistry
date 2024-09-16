#ifndef MOLHEADERDEF
#define MOLHEADERDEF

#include "../libgf/gf.hpp"

#include <fstream>
#include <iostream>
#include <string>

class Molecule{
	public:
		Molecule(std::string file, std::string bfs);

		int Natoms;
		std::vector<int> Zvals;
		int charge;
		int Nelec;
		int NUPDOWN;
		std::string basis;
		std::vector<std::vector<double>> xyz;
		std::vector<GF> AOs;
};

std::vector<GF> AOfunctions(std::string bfs, int Zval, std::vector<double> xyz);

#endif
