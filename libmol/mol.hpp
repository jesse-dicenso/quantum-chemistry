#ifndef MOLHEADERDEF
#define MOLHEADERDEF

#include "../libgf/gf.hpp"

#include <fstream>
#include <iostream>
#include <string>

class Molecule{
	public:
		Molecule(const std::string& file, const std::string& bfs);

		int Natoms;
		std::vector<int> Zvals;
		bool heteronuclear;
		int charge;
		int Nelec;
		int NUPDOWN;
		std::string basis;
		std::vector<std::vector<double>> xyz;
		std::vector<GF> AOs;
};

std::vector<GF> AOfunctions(const std::string& bfs, int Zval, const std::vector<double>& xyz, int atom_idx);
extern const std::vector<std::string> elements;

#endif
