#include "mol.hpp"

Molecule::Molecule(std::string file, std::string bfs){
	std::ifstream inpfile(file);
	assert(inpfile.good());

	Nelec = 0;
	inpfile >> Natoms >> charge >> NUPDOWN;
	Nelec = -charge;
	heteronuclear = false;
	Zvals.resize(Natoms);
	xyz.resize(Natoms);
	for(int i = 0; i < Natoms; i++){
		xyz[i].resize(3);
		inpfile >> Zvals[i] >> xyz[i][0] >> xyz[i][1] >> xyz[i][2];
		if(Zvals[i] != Zvals[0]){
			heteronuclear = true;
		}
		Nelec += Zvals[i];
	}
	inpfile.close();

	basis = bfs;
	std::vector<GF> temp;	
	for(int j = 0; j < Natoms; j++){
		temp = AOfunctions(basis, Zvals[j], xyz[j], j);
		for(int k = 0; k < temp.size(); k++){
			AOs.push_back(temp[k]);
		}
	}
}

std::vector<GF> AOfunctions(std::string bfs, int Zval, const std::vector<double>& pos, int atom_idx){
	// Element symbols ( elements[i] == elements[Zval-1] )
	std::vector<std::string> elements = {
    		"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    		"Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    		"Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    		"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    		"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    		"Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    		"Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    		"Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    		"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    		"Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    		"Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    		"Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
	};
	std::string element_symbol = elements[Zval-1];

	// Cartesian angular momenta up to l=2; ( Ncl = (l+1)(l+2)/2 )
	std::vector<int> Ls, Lpx, Lpy, Lpz, Ldx2, Ldy2, Ldz2, Ldxy, Ldyz, Ldzx;
	Ls   = {0,0,0};

	Lpx  = {1,0,0}; 
	Lpy  = {0,1,0}; 
	Lpz  = {0,0,1};
 
	Ldx2 = {2,0,0};
	Ldy2 = {0,2,0};
	Ldz2 = {0,0,2};
	Ldxy = {1,1,0};
	Ldyz = {0,1,1};
	Ldzx = {1,0,1};
	
	// Read in basis set and construct AO functions
	assert((bfs=="STO-3G") || (bfs=="def2-SVP"));
	std::ifstream bfsfile("../libmol/"+bfs+"/"+element_symbol);
	assert(bfsfile.good());
	
	std::vector<GF> orbitals;

	int numshells = 0;
	bfsfile >> numshells;
	for(int i = 0; i < numshells; i++){
		std::string shell;
		int clen;
		double zeta_scale;

		bfsfile >> shell >> clen >> zeta_scale;
		assert((shell=="S") || (shell=="SP") || (shell=="P") || (shell=="D") || (shell=="F"));

		std::vector<double> zeta(clen), d(clen);

		if(shell=="S"){
			for(int j = 0; j < clen; j++){
				bfsfile >> zeta[j] >> d[j];
				zeta[j]*=zeta_scale;
			}
			GF S(zeta, d, pos, Ls, atom_idx);
			orbitals.push_back(S);
		}
		else if(shell=="SP"){
			std::vector<double> d2(clen);
			for(int j = 0; j < clen; j++){
				bfsfile >> zeta[j] >> d[j] >> d2[j];
				zeta[j]*=zeta_scale;
			}
			GF S (zeta, d , pos, Ls , atom_idx);
			GF Px(zeta, d2, pos, Lpx, atom_idx);
			GF Py(zeta, d2, pos, Lpy, atom_idx);
			GF Pz(zeta, d2, pos, Lpz, atom_idx);
			orbitals.push_back(S );	
			orbitals.push_back(Px);	
			orbitals.push_back(Py);	
			orbitals.push_back(Pz);	
		}
		else if(shell=="P"){
			for(int j = 0; j < clen; j++){
				bfsfile >> zeta[j] >> d[j];
				zeta[j]*=zeta_scale;
			}
			GF Px(zeta, d, pos, Lpx, atom_idx);
			GF Py(zeta, d, pos, Lpy, atom_idx);
			GF Pz(zeta, d, pos, Lpz, atom_idx);
			orbitals.push_back(Px);	
			orbitals.push_back(Py);	
			orbitals.push_back(Pz);	
		}
		else if(shell=="D"){
			for(int j = 0; j < clen; j++){
				bfsfile >> zeta[j] >> d[j];
				zeta[j]*=zeta_scale;
			}
			GF Dx2(zeta, d, pos, Ldx2, atom_idx);
			GF Dy2(zeta, d, pos, Ldy2, atom_idx);
			GF Dz2(zeta, d, pos, Ldz2, atom_idx);
			GF Dxy(zeta, d, pos, Ldxy, atom_idx);
			GF Dyz(zeta, d, pos, Ldyz, atom_idx);
			GF Dzx(zeta, d, pos, Ldzx, atom_idx);
			orbitals.push_back(Dx2);	
			orbitals.push_back(Dy2);	
			orbitals.push_back(Dz2);	
			orbitals.push_back(Dxy);	
			orbitals.push_back(Dyz);	
			orbitals.push_back(Dzx);	
		}
	}	

	bfsfile.close();
	
	return orbitals;
}
