#include "hartree-fock.hpp"

using namespace std;

Matrix density_matrix(Matrix C, int N);
double E0(Matrix P, Matrix Hcore, Matrix F);

int main(){
	cout << fixed;
	cout << setprecision(10);
	cout << scientific;

	string input;
	const double eps = 1e-6;
	const int max_cycles = 50;

	cout << "Input file:\n";
	cin >> input;
	//input = "H2O.inp";

	int N;
	int K;
	int cycles = 1;
	double err = eps + 1;
	vector<Matrix> temp_e_c;
	vector<vector<vector<vector<double>>>> eris;

	double nuc;
	double Eo;
	double temp_Eo;
	double E_tot;

	//////////////////////////////////////////////////
	// s		// overlap			//
	// t		// kinetic			//
	// v		// nuclear attraction		//
	// hcore	// core Hamiltonian		//
	// g		// two-electron fock component	//
	// x		// orthogonalization (symmetric)//
	// p		// density			//
	// c		// coefficient			//
	// c		// orthogonalized coefficient	//
	// f		// fock				//
	// fo		// orthogonalized fock		//
	// e		// orbital energies		//
	//////////////////////////////////////////////////
	
	cout << "********************\n";
	cout << "*                  *\n";
	cout << "*   HARTREE-FOCK   *\n";
	cout << "*                  *\n";
	cout << "********************\n\n";
	cout << "JESSE DICENSO   \n";
	cout << "SUMMER 2024     \n\n";

	cout << "Reading input file...\n\n";
	
	Molecule M(input);
	N = M.Nelec;
	K = M.AOs.size();

	cout << "basis      STO-3G\n";
	cout << "eps        " << eps << '\n';
	cout << "max_cycles " << max_cycles << "\n\n";
	cout << "--------------------------------------------------------------\n";
	cout << "Z" << setw(20) << "x (Bohr)" << setw(20) << "y (Bohr)" << setw(20) <<  "z (Bohr)" << '\n';
	cout << "--------------------------------------------------------------\n";
	for(int i = 0; i < M.Zvals.size(); i++){
		cout << M.Zvals[i] << setw(20) << M.xyz[i][0] << setw(20) << M.xyz[i][1] << setw(20) << M.xyz[i][2] << '\n';
	}
	cout << "--------------------------------------------------------------\n";
	cout << "Success! There are " << N << " electrons and " << K << " basis functions.\n\n"; 
	
	nuc = nucrepl(M.Zvals, M.xyz);	

	cout << "Nuclear Repulsion Energy = " << nuc << " Ha\n\n";
	
	eris = ERIs(M.AOs);

	Matrix s = overlap(M.AOs);
	
	cout << "OVERLAP\n";
	s.printMatrix();
	
	Matrix t = kinetic(M.AOs);
	
	cout << "KINETIC ENERGY\n";
	t.printMatrix();
	
	Matrix v = nuclear(M.AOs, M.Zvals, M.xyz);
	
	cout << "NUCLEAR ATTRACTION\n";
	v.printMatrix();
	
	Matrix hcore = t + v;
	
	cout << "CORE HAMILTONIAN\n";
	hcore.printMatrix();

	cout << "Using core Hamiltonian as initial guess to Fock matrix.\n\n";
	cout << "------------------------------------------\n\n";
	cout << "              *************\n";
	cout << "              * BEGIN SCF *\n";
	cout << "              *************\n\n";
	cout << "------------------------------------------\n";
	cout << "N" << setw(20) << "E" <<  setw(20) << "err" << '\n';
	cout << "------------------------------------------\n";

	Matrix x = inv_sqrt(s);

	Matrix p = zero(K,K);
	Matrix g = G(p, eris);
	Matrix f = hcore + g;
	
	Matrix fo = transpose(x) * f * x;
	temp_e_c = QR_diagonalize(fo);

	Matrix e = temp_e_c[0];
	Matrix co = temp_e_c[1];

	Matrix c = x * co;
	p = density_matrix(c, N);
	
	while((abs(err) > eps) && (cycles <= max_cycles)){
		g = G(p, eris);
		f = hcore + g;
		temp_Eo = E0(p, hcore, f);
		err = temp_Eo - Eo;
		Eo = temp_Eo;

		fo = transpose(x) * f * x;
		temp_e_c = QR_diagonalize(fo);
		e = temp_e_c[0];
		co = temp_e_c[1];
		
		c = x * co;
		p = density_matrix(c, N);
		
		cout << cycles << setw(20) << Eo << setw(20) << err << '\n';
		cycles+=1;
	}
	cout << "------------------------------------------\n";
	cout << '\n';

	if(cycles > max_cycles){
		cout << "No convergence after " << max_cycles << " cycles!" << '\n';
	}
	else{
		E_tot = Eo + nuc;

		cout << "Total energy (Ha) = " << Eo + nuc << "\n\n";

		cout << "MO coefficients\n";
		c.printMatrix();
		
		cout << "Orbital Energies\n";
		e.printMatrix();
	}
}

Matrix density_matrix(Matrix C, int N){
	Matrix P(C.rows, C.cols);
	for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			double sum = 0;
			for(int a = 0; a < N/2; a++){
				sum += C.matrix[i][a]*C.matrix[j][a];
			} 
			P.matrix[i][j] = 2*sum;
		}
	}
	return P;
}

double E0(Matrix P, Matrix Hcore, Matrix F){
	double sum = 0;
	for(int i = 0; i < P.rows; i++){
		for(int j = 0; j < P.cols; j++){
			sum += P.matrix[j][i]*(Hcore.matrix[i][j]+F.matrix[i][j]);
		}
	}
	return 0.5 * sum;
}
