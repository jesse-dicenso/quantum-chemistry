#include "hartree-fock.hpp"

using namespace std;

vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S);
vector<double> Mulliken_PA(Molecule M, Matrix P, Matrix S);
void print_orbitals(Matrix E, Matrix C, int N, int K);

int main(int argc, char* argv[]){
	cout << fixed;
	cout << scientific;

	ofstream out("outfile.dat");
	auto coutbuf = cout.rdbuf(out.rdbuf());

	string infile;
	int sps;
	double eps;
	int max_cycles;
	string pop;

	cin >> infile;
	cin >> sps;
	cin >> eps;
	cin >> max_cycles;
	cin >> pop;

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

	cout << "Reading input...\n\n";
	
	Molecule M(infile);
	N = M.Nelec;
	K = M.AOs.size();
	cout << "\tinfile\t\t" << infile << '\n';
	cout << "\tbasis\t\tSTO-3G\n";
	cout << "\tsps\t\t" << sps << ' ';
	if(sps==0){
		cout << "(fixed-point)\n";
	}
	else if(sps>0){
		cout << "(DIIS)\n";
	}
	else{
		cerr << "ERR: Unknown SCF algorithm\n";
		return 0;
	}
	cout << "\teps\t\t" << setprecision(2) << eps << '\n';
	cout << "\tmax_cycles\t" << max_cycles << '\n';
	cout << "\tpopulation\t" << pop << "\n\n";
	cout << setprecision(10);
	cout << "--------------------------------------------------------------\n";
	cout << "Z" << setw(20) << "x (Bohr)" << setw(20) << "y (Bohr)" << setw(20) <<  "z (Bohr)" << '\n';
	cout << "--------------------------------------------------------------\n";
	for(int i = 0; i < M.Zvals.size(); i++){
		cout << M.Zvals[i] << setw(20) << M.xyz[i][0] << setw(20) << M.xyz[i][1] << setw(20) << M.xyz[i][2] << '\n';
	}
	cout << "--------------------------------------------------------------\n";
	cout << "Success! There are " << N << " electrons and " << K << " basis functions.\n\n"; 
	
	nuc = nucrepl(M.Zvals, M.xyz);	

	cout << "Nuclear Repulsion Energy (Ha) = " << nuc << '\n';
	
	eris = ERIs(M.AOs);
	
	Matrix s = overlap(M.AOs);
	Matrix hcore = kinetic(M.AOs) + nuclear(M.AOs, M.Zvals, M.xyz);

	cout << "Using core Hamiltonian as initial guess to Fock matrix.\n\n";
	cout << "              *************\n";
	cout << "              * BEGIN SCF *\n";
	cout << "              *************\n\n";
	cout << "----------------------------------------------------\n";
	cout << setw(3) << "N" << setw(20) << "E" <<  setw(20) << "err" << setw(10) << "step\n";
	cout << "----------------------------------------------------\n";

	Matrix x = m_inv_sqrt(s);
	Matrix p(K, K);
	Matrix g(K, K);
	Matrix f = hcore;
	Eo = E0(p, hcore, f);
	Matrix fo = transpose(x) * f * x;
	temp_e_c = diagonalize(fo);
	Matrix e = temp_e_c[0];
	Matrix co = temp_e_c[1];
	Matrix c = x * co;
	p = density_matrix(c, N);

	cout.flush();	
	
	if(sps==0){
		while((abs(err) > eps) && (cycles <= max_cycles)){
			FPI(hcore, eris, x, &p, &g, &f, &fo, &e, &co, &c, &Eo, &err, N);
			cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err << setw(10) << "fp\n";
			cout.flush();
			cycles+=1;
		}
	}
	else{
		vector<Matrix> SPf(sps, zero(K, K));
		vector<Matrix> SPe(sps, zero(K, K));
		while((abs(err) > eps) && (cycles <= max_cycles)){
			DIIS(s, hcore, eris, x, &p, &g, &f, &fo, &e, &co, &c, &Eo, &err, N, cycles, SPf.data(), SPe.data(), sps);
			cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err;
			if(cycles < sps){
				cout << setw(10) << "fp\n";
			}
			else{
				cout << setw(10) << "diis\n";
			}
			cout.flush();
			cycles+=1;
		}

	}
	cout << "----------------------------------------------------\n";

	if(cycles > max_cycles){
		cerr << "ERR: No convergence after " << max_cycles << " cycles!\n";
		return 0;
	}
	else{
		cout << "Convergence criterion met; exiting SCF loop.\n\n";
		E_tot = Eo + nuc;

		cout << "Total energy (Ha) = " << Eo + nuc << "\n\n";
		
		print_orbitals(e, c, N, K);

		vector<double> pa(M.Zvals.size());
		if(pop=="lowdin"){
			pa = Lowdin_PA(M, p, s);
			cout << "=======================\n";
			cout << "Löwdin Pop. Analysis\n";
		}
		else if(pop=="mulliken"){
			pa = Mulliken_PA(M, p, s);
			cout << "=======================\n";
			cout << "Mulliken Pop. Analysis\n";
		}
		cout << "=======================\n";
		cout << setw(3) << "idx" << setw(21) << "charge\n";
		cout << "=======================\n";
		double sum_chg = 0;
		cout << fixed;
		cout << setprecision(10);
		for(int i = 0; i < M.Zvals.size(); i++){
			cout << setw(3) << i+1 << setw(20) << pa[i] << '\n';
			sum_chg += pa[i];
		}
		cout << "=======================\n";
		cout << "Sum of atomic charges = " << sum_chg << "\n\n";

		cout << "\n***************************************\n";
		cout <<   "*                                     *\n";
		cout <<   "*     Thank you, have a nice day!     *\n";
		cout <<   "*                                     *\n";
		cout <<   "***************************************\n\n";

	}
}

vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S){
	vector<double> LPA(M.Zvals.size());
	Matrix ShPSh = m_sqrt(S) * P * m_sqrt(S);
	for(int i = 0; i < M.Zvals.size(); i++){
		double sum = 0;
		vector<int> positions;
		for(int j = 0; j < M.AOs.size(); j++){
			if(M.AOs[j].xyz==M.xyz[i]){
				sum += ShPSh.matrix[j][j];
			}
		}
		LPA[i] = M.Zvals[i] - sum;
	}
	return LPA;
}

vector<double> Mulliken_PA(Molecule M, Matrix P, Matrix S){
	vector<double> MPA(M.Zvals.size());
	Matrix PS = P * S;
	for(int i = 0; i < M.Zvals.size(); i++){
		double sum = 0;
		vector<int> positions;
		for(int j = 0; j < M.AOs.size(); j++){
			if(M.AOs[j].xyz==M.xyz[i]){
				sum += PS.matrix[j][j];
			}
		}
		MPA[i] = M.Zvals[i] - sum;
	}
	return MPA;
}

void print_orbitals(Matrix E, Matrix C, int N, int K){
	cout << "********************\n";
	cout << "* Orbital Energies *\n";
	cout << "********************\n\n";
	cout << "Occupied:\n";
	for(int i = 0; i < N/2; i++){
		cout << setw(4) << 'E' << i+1 << setw(20) << E.matrix[i][i] << '\n';
	}
	if((K-N/2)>0){
		cout << "Virtual:\n";
		for(int i = N/2; i < K; i++){
			cout << setw(4) << 'E' << i+1 << setw(20) << E.matrix[i][i] << '\n'; 
		}
	}
	cout << '\n';
	cout << "************************\n";
	cout << "* Orbital Coefficients *\n";
	cout << "************************\n\n";
	cout << "Occupied:\n";
	int count = 0;
	for(int i = 0; i < N/2; i++){
		cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
		for(int j = 0; j < K; j++){
			cout << C.matrix[j][i] << setw(20);
			count++;
			if(count>=5){
				cout << '\n' << setw(25);
				count = 0;
			}
		}
		count = 0;
		cout << '\n';
	}
	if((K-N/2)>0){
		cout << "Virtual:\n";
		for(int i = N/2; i < K; i++){
			cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
			for(int j = 0; j < K; j++){
				cout << C.matrix[j][i] << setw(20);
				count++;
				if(count>=5){
					cout << '\n' << setw(25);
					count = 0;
			}
		}
		count = 0;
		cout << '\n';
		}
	}
	cout << '\n';
}
