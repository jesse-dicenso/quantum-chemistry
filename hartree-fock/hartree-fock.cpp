#include "hartree-fock.hpp"

using namespace std;

Matrix density_matrix(Matrix C, int N){
	Matrix P;	
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

int main(){
	cout << fixed;
	cout << setprecision(15);
	const double eps = 1e-6;
	const int max_cycles = 50;
	double err;

	int N;
	int K;
	int cycles;
	vector<Matrix> temp_e_c;
	vector<vector<vector<vector<double>>>> eris;

	double nuc;
	double Eo;
	double temp_Eo;

	//////////////////////////////////////////////////
	Matrix s;	// overlap			//
	Matrix t;	// kinetic			//
	Matrix v;	// nuclear attraction		//
	Matrix hcore;	// core Hamiltonian		//
	Matrix g;	// two-electron fock component	//
	Matrix x;	// orthogonalization		//
	Matrix p;	// density			//
	Matrix c;	// coefficient			//
	Matrix co;	// orthogonalized coefficient	//
	Matrix f;	// fock				//
	Matrix fo;	// orthogonalized fock		//
	Matrix e;	// orbital energies		//
							//
	//////////////////////////////////////////////////
	
	cout << "********************\n";
	cout << "*                  *\n";
	cout << "*   HARTREE-FOCK   *\n";
	cout << "*                  *\n";
	cout << "********************\n\n";
	cout << "JESSE DICENSO   \n";
	cout << "SUMMER 2024     \n\n";

	cout << "Reading input file...\n\n";
	
	Molecule M("h2.inp");
	N = M.Nelec;
	K = M.AOs.size();

	cout << "basis      STO-3G\n";
	cout << "eps        " << eps << '\n';
	cout << "max_cycles " << max_cycles << "\n\n";

	cout << "-------------------------------------------------------------\n";
	cout << "Z" << setw(20) << "x (Bohr)" << setw(20) << "y (Bohr)" << setw(20) << "z (Bohr)" << '\n';
	cout << "-------------------------------------------------------------\n";
	for(int i = 0; i < M.Zvals.size(); i++){
		cout << M.Zvals[i] << setw(20) << M.xyz[i][0] << setw(20) << M.xyz[i][1] << setw(20) << M.xyz[i][2] << '\n';
	}
	cout << "-------------------------------------------------------------\n";
	cout << "Success! There are " << N << " electrons and " << K << " basis functions.\n\n"; 
	
	nuc = nucrepl(M.Zvals, M.xyz);	

	cout << "Nuclear Repulsion Energy = " << nuc << " Ha\n\n";
	
	s = overlap(M.AOs);
	cout << "OVERLAP\n";
	s.printMatrix();
	cout << '\n';

	t = kinetic(M.AOs);
	cout << "KINETIC ENERGY\n";
	t.printMatrix();
	cout << '\n';
	
	v = nuclear(M.AOs, M.Zvals, M.xyz);
	cout << "NUCLEAR ATTRACTION\n";
	v.printMatrix();
	cout << '\n';

	hcore = t + v;
	cout << "CORE HAMILTONIAN\n";
	hcore.printMatrix();
	cout << '\n';
	
	eris = ERIs(M.AOs);

	x = inv_sqrt(s);

	// Use Hcore as initial guess for P
	p = hcore;

	// SCF cycle
	g = G(p, eris);
	f = hcore + g;
	fo = transpose(x) * f * x;

	cout << "CORE HAMILTONIAN\n";
	fo.printMatrix();
	cout << '\n';
	
	temp_e_c = QR_diagonalize(fo);

	cout << "diagged\n";

	e = temp_e_c[0];
	co = temp_e_c[1];
	c = x * co;

	p = density_matrix(c, N);		
		
	Eo = E0(p, hcore, f);
	cycles = 1;
	while((abs(err) > eps) && (cycles <= max_cycles)){
		g = G(p, eris);
		f = hcore + g;
		fo = transpose(x) * f * x;
		
		temp_e_c = QR_diagonalize(fo);
		e = temp_e_c[0];
		co = temp_e_c[1];
		c = x * co;

		p = density_matrix(c, N);		
		
		temp_Eo = E0(p, hcore, f);
		err = temp_Eo - Eo;
		Eo = temp_Eo;
		cycles+=1;
	}
	cout << "Total energy (Ha) = " << Eo + nuc << '\n';
}
