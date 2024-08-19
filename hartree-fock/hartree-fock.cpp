#include "hartree-fock.hpp"
#include <bits/stdc++.h>

using namespace std;

Matrix density_matrix(Matrix C, int N);
double E0(Matrix P, Matrix Hcore, Matrix F);
vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S);
vector<double> Mulliken_PA(Molecule M, Matrix P, Matrix S);

vector<int> e_order(Matrix E);
Matrix e_reorder(Matrix E, vector<int> idx);
Matrix c_reorder(Matrix C, vector<int> idx);

int main(int argc, char* argv[]){
	cout << fixed;
	cout << scientific;

	ofstream out("outfile.dat");
	auto coutbuf = cout.rdbuf(out.rdbuf());

	string infile;
	double eps;
	int max_cycles;
	string pop;

	cin >> infile;
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
	cout << "SUMMER 2024     \n\n\n";

	cout << "Reading input file...\n\n";
	
	Molecule M(infile);
	N = M.Nelec;
	K = M.AOs.size();

	cout << "\tbasis      STO-3G\n";
	cout << "\teps        " << setprecision(2) << eps << '\n';
	cout << "\tmax_cycles " << max_cycles << '\n';
	cout << "\tpopulation " << pop << "\n\n";
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
	cout << "              *************\n";
	cout << "              * BEGIN SCF *\n";
	cout << "              *************\n\n";
	cout << "--------------------------------------------\n";
	cout << setw(3) << "N" << setw(20) << "E" <<  setw(20) << "err" << '\n';
	cout << "--------------------------------------------\n";

	Matrix x = inv_sqrt(s);

	Matrix p = zero(K,K);
	Matrix g = zero(K,K);
	Matrix f = hcore;
	Eo = E0(p, hcore, f);
	
	Matrix fo = transpose(x) * f * x;
	temp_e_c = QR_diagonalize(fo);

	Matrix e = temp_e_c[0];
	vector e_o = e_order(e);
	e = e_reorder(e, e_o);

	Matrix co = temp_e_c[1];

	Matrix c = x * co;
	c = c_reorder(c, e_o);

	p = density_matrix(c, N);

	cout.flush();	
	while((abs(err) > eps) && (cycles <= max_cycles)){
		g = G(p, eris);
		f = hcore + g;
		temp_Eo = E0(p, hcore, f);
		err = temp_Eo - Eo;
		Eo = temp_Eo;

		fo = transpose(x) * f * x;
		temp_e_c = QR_diagonalize(fo);
		e = temp_e_c[0];
		e_o = e_order(e);
		e = e_reorder(e, e_o);
		co = temp_e_c[1];
		
		c = x * co;
		c = c_reorder(c, e_o);
		p = density_matrix(c, N);
	
		cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err << '\n';
		cout.flush();
		cycles+=1;
	}
	cout << "------------------------------------------\n";

	if(cycles > max_cycles){
		cout << "No convergence after " << max_cycles << " cycles!" << '\n';
	}
	else{
		cout << "Convergence criterion met; exiting SCF loop.\n\n";
		E_tot = Eo + nuc;

		cout << "Total energy (Ha) = " << Eo + nuc << "\n\n";

		cout << "MO coefficients\n";
		c.printMatrix();
		
		cout << "Orbital Energies\n";
		e.printMatrix();

		cout << "Final Density Matrix\n";
		p.printMatrix();

		vector<double> pa(M.Zvals.size());
		if(pop=="lowdin"){
			pa = Lowdin_PA(M, p, s);
			cout << "Löwdin Population Analysis\n";
		}
		else if(pop=="mulliken"){
			pa = Mulliken_PA(M, p, s);
			cout << "Mulliken Population Analysis\n";
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

vector<int> e_order(Matrix E){
	vector<int> list_idx(E.cols);
	vector<double> list_e(E.cols);
	for(int i = 0; i < E.cols; i++){
		list_e[i] = E.matrix[i][i];
	}
	sort(list_e.begin(), list_e.end());
	for(int j = 0; j < E.cols; j++){
		for(int k = 0; k < E.cols; k++){
			if(list_e[j] == E.matrix[k][k]){
				list_idx[j] = k;
			}
		}
	}
	return list_idx;
}

Matrix e_reorder(Matrix E, vector<int> idx){
	Matrix e_ordered(E.rows, E.cols);
	for(int i = 0; i < E.cols; i++){
		e_ordered.matrix[i][i] = E.matrix[idx[i]][idx[i]];
	}
	return e_ordered;
}

Matrix c_reorder(Matrix C, vector<int> idx){
	Matrix c_ordered(C.rows, C.cols);
	for(int i = 0; i < C.cols; i++){
		for(int j = 0; j < C.rows; j++){
			c_ordered.matrix[j][i] = C.matrix[j][idx[i]];
		}
	}
	return c_ordered;
}

vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S){
	vector<double> LPA(M.Zvals.size());
	Matrix ShPSh = sqrt(S) * P * sqrt(S);
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
