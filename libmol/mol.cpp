#include "mol.hpp"

Molecule::Molecule(std::string file, std::string bfs){
	std::ifstream read(file);
	assert(read.good());

	Nelec = 0;
	read >> Natoms >> charge >> NUPDOWN;
	Nelec = -charge;
	Zvals.resize(Natoms);
	xyz.resize(Natoms);
	for(int i = 0; i < Natoms; i++){
		xyz[i].resize(3);
		read >> Zvals[i] >> xyz[i][0] >> xyz[i][1] >> xyz[i][2];
		Nelec += Zvals[i];
	}
	read.close();

	basis = bfs;
	std::vector<GF> temp;	
	for(int j = 0; j < Natoms; j++){
		temp = AOfunctions(basis, Zvals[j], xyz[j]);
		for(int k = 0; k < temp.size(); k++){
			AOs.push_back(temp[k]);
		}
	}
}

std::vector<GF> AOfunctions(std::string bfs, int Zval, std::vector<double> pos){
	std::vector<int> Ls ({0,0,0});
	std::vector<int> Lpx({1,0,0});
	std::vector<int> Lpy({0,1,0});
	std::vector<int> Lpz({0,0,1});
	// only STO-3G available now...sorry
	assert(bfs=="STO-3G");
	std::vector<GF> orbitals;
		// * H * //
		if(Zval==1){
			// 1s
			std::vector<double> a1s({0.3425250914e1,0.6239137298,0.1688554040});
			std::vector<double> d1s({0.1543289673,0.5353281423,0.4446345422});
			GF H1s(a1s, d1s, pos, Ls);
			orbitals.push_back(H1s);
		}
		// * He * //
		else if(Zval==2){
			// 1s
			// Basis set exchange:
			//std::vector<double> a1s({0.6362421394e1,0.1158922999e1,0.3136497915});
			//std::vector<double> d1s({0.1543289673,0.5353281423,0.4446345422});
			// Szabo-Ostlund:
			std::vector<double> a1s({0.480844,1.776691,9.753934});
			std::vector<double> d1s({0.444635,0.535328,0.154329});
			GF He1s(a1s, d1s, pos, Ls);
			orbitals.push_back(He1s);
		}
		// * C * //
		else if(Zval==6){
			// 1s
			std::vector<double> a1s({0.7161683735e2,0.1304509632e2,0.3530512160e1});
			std::vector<double> d1s({0.1543289673,0.5353281423,0.4446345422});
			GF C1s(a1s, d1s, pos, Ls);
			orbitals.push_back(C1s);
			// 2s, 2p
			std::vector<double> a2sp({0.2941249355e1,0.6834830964,0.2222899159});
			std::vector<double> d2s ({-0.9996722919e-1,0.3995128261,0.7001154689});
			std::vector<double> d2p ({0.1559162750,0.6076837186,0.3919573931});
			GF C2s (a2sp, d2s, pos, Ls );
			GF C2px(a2sp, d2p, pos, Lpx);
			GF C2py(a2sp, d2p, pos, Lpy);
			GF C2pz(a2sp, d2p, pos, Lpz);
			orbitals.push_back(C2s );
			orbitals.push_back(C2px);
			orbitals.push_back(C2py);
			orbitals.push_back(C2pz);
		}
		// * O * //
		else if(Zval==8){
			// 1s
			std::vector<double> a1s({0.1307093214e3,0.2380886605e2,0.6443608313e1});
			std::vector<double> d1s({0.1543289673,0.5353281423,0.4446345422});
			GF O1s(a1s, d1s, pos, Ls);
			orbitals.push_back(O1s);
			// 2s, 2p
			std::vector<double> a2sp({0.5033151319e1,0.1169596125e1,0.3803889600});
			std::vector<double> d2s ({-0.9996722919e-1,0.3995128261,0.7001154689});
			std::vector<double> d2p ({0.1559162750,0.6076837186,0.3919573931});
			GF O2s (a2sp, d2s, pos, Ls );
			GF O2px(a2sp, d2p, pos, Lpx);
			GF O2py(a2sp, d2p, pos, Lpy);
			GF O2pz(a2sp, d2p, pos, Lpz);
			orbitals.push_back(O2s );
			orbitals.push_back(O2px);
			orbitals.push_back(O2py);
			orbitals.push_back(O2pz);
		}
	return orbitals;
}
