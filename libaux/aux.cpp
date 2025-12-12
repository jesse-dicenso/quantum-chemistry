#include "aux.hpp"

std::vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S){
	std::vector<double> LPA(M.Zvals.size());
	Matrix ShPSh = m_sqrt(S) * P * m_sqrt(S);
	for(int i = 0; i < M.Zvals.size(); i++){
		double sum = 0;
		for(int j = 0; j < M.AOs.size(); j++){
			if(M.AOs[j].atom_index==i){
				sum += ShPSh.matrix[j][j];
			}
		}
		LPA[i] = M.Zvals[i] - sum;
	}
	return LPA;
}

std::vector<double> Mulliken_PA(Molecule M, Matrix P, Matrix S){
	std::vector<double> MPA(M.Zvals.size());
	Matrix PS = P * S;
	for(int i = 0; i < M.Zvals.size(); i++){
		double sum = 0;
		for(int j = 0; j < M.AOs.size(); j++){
			if(M.AOs[j].atom_index==i){
				sum += PS.matrix[j][j];
			}
		}
		MPA[i] = M.Zvals[i] - sum;
	}
	return MPA;
}

void R_print_orbitals(Matrix E, Matrix C, int Nocc, int Kb){
	std::cout << "***************\n";
	std::cout << "* MO Energies *\n";
	std::cout << "***************\n\n";
	std::cout << "Occupied:\n";
	for(int i = 0; i < Nocc/2; i++){
		std::cout << std::setw(4) << 'E' << i+1 << std::setw(20) << E.matrix[i][i] << std::endl;
	}
	if(Kb>(Nocc/2)){
		std::cout << "Virtual:\n";
		for(int i = Nocc/2; i < Kb; i++){
			std::cout << std::setw(4) << 'E' << i+1 << std::setw(20) << E.matrix[i][i] << std::endl; 
		}
	}
	std::cout << '\n';
	std::cout << "*******************\n";
	std::cout << "* MO Coefficients *\n";
	std::cout << "*******************\n\n";
	std::cout << "Occupied:\n";
	int count = 0;
	for(int i = 0; i < Nocc/2; i++){
		std::cout << std::setw(5) << "MO" << i+1 << '\n' << std::setw(25); 
		for(int j = 0; j < Kb; j++){
			std::cout << C.matrix[j][i] << std::setw(20);
			count++;
			if(count>=5){
				std::cout << '\n' << std::setw(25);
				count = 0;
			}
		}
		count = 0;
		std::cout << '\n';
	}
	if(Kb>(Nocc/2)){
		std::cout << "Virtual:\n";
		for(int i = Nocc/2; i < Kb; i++){
			std::cout << std::setw(5) << "MO" << i+1 << '\n' << std::setw(25); 
			for(int j = 0; j < Kb; j++){
				std::cout << C.matrix[j][i] << std::setw(20);
				count++;
				if(count>=5){
					std::cout << '\n' << std::setw(25);
					count = 0;
			}
		}
		count = 0;
		std::cout << '\n';
		}
	}
	std::cout << '\n';
}

void UR_print_orbitals(Matrix Ea, Matrix Eb, Matrix Ca, Matrix Cb, int Nocca, int Noccb, int Kb){
	std::cout << "***************\n";
	std::cout << "* MO Energies *\n";
	std::cout << "***************\n\n";
	std::cout << "Occupied (alpha):\n";
	for(int i = 0; i < Nocca; i++){
		std::cout << std::setw(4) << 'E' << i+1 << std::setw(20) << Ea.matrix[i][i] << '\n';
	}
	if(Kb>Nocca){
		std::cout << "Virtual (alpha):\n";
		for(int i = Nocca; i < Kb; i++){
			std::cout << std::setw(4) << 'E' << i+1 << std::setw(20) << Ea.matrix[i][i] << '\n'; 
		}
	}
	std::cout << '\n';
	std::cout << "Occupied (beta):\n";
	for(int i = 0; i < Noccb; i++){
		std::cout << std::setw(4) << 'E' << i+1 << std::setw(20) << Eb.matrix[i][i] << '\n';
	}
	if(Kb>Noccb){
		std::cout << "Virtual (beta):\n";
		for(int i = Noccb; i < Kb; i++){
			std::cout << std::setw(4) << 'E' << i+1 << std::setw(20) << Eb.matrix[i][i] << '\n'; 
		}
	}
	std::cout << '\n';
	std::cout << "*******************\n";
	std::cout << "* MO Coefficients *\n";
	std::cout << "*******************\n\n";
	std::cout << "Occupied (alpha):\n";
	int count = 0;
	for(int i = 0; i < Nocca; i++){
		std::cout << std::setw(5) << "MO" << i+1 << '\n' << std::setw(25); 
		for(int j = 0; j < Kb; j++){
			std::cout << Ca.matrix[j][i] << std::setw(20);
			count++;
			if(count>=5){
				std::cout << '\n' << std::setw(25);
				count = 0;
			}
		}
		count = 0;
		std::cout << '\n';
	}
	if(Kb>Nocca){
		std::cout << "Virtual (alpha):\n";
		for(int i = Nocca; i < Kb; i++){
			std::cout << std::setw(5) << "MO" << i+1 << '\n' << std::setw(25); 
			for(int j = 0; j < Kb; j++){
				std::cout << Ca.matrix[j][i] << std::setw(20);
				count++;
				if(count>=5){
					std::cout << '\n' << std::setw(25);
					count = 0;
			}
		}
		count = 0;
		std::cout << '\n';
		}
	}
	std::cout << '\n';
	std::cout << "Occupied (beta):\n";
	count = 0;
	for(int i = 0; i < Noccb; i++){
		std::cout << std::setw(5) << "MO" << i+1 << '\n' << std::setw(25); 
		for(int j = 0; j < Kb; j++){
			std::cout << Cb.matrix[j][i] << std::setw(20);
			count++;
			if(count>=5){
				std::cout << '\n' << std::setw(25);
				count = 0;
			}
		}
		count = 0;
		std::cout << '\n';
	}
	if(Kb>Noccb){
		std::cout << "Virtual (beta):\n";
		for(int i = Noccb; i < Kb; i++){
			std::cout << std::setw(5) << "MO" << i+1 << '\n' << std::setw(25); 
			for(int j = 0; j < Kb; j++){
				std::cout << Cb.matrix[j][i] << std::setw(20);
				count++;
				if(count>=5){
					std::cout << '\n' << std::setw(25);
					count = 0;
			}
		}
		count = 0;
		std::cout << '\n';
		}
	}
	std::cout << '\n';
}
