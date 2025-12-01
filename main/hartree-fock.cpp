#include "hartree-fock.hpp"

using namespace std;

vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S);
vector<double> Mulliken_PA(Molecule M, Matrix P, Matrix S);
void R_print_orbitals(Matrix E, Matrix C, int Nocc, int Kb);
void UR_print_orbitals(Matrix Ea, Matrix Eb, Matrix Ca, Matrix Cb, int Nocca, int Noccb, int Kb);

int main(int argc, char* argv[]){
	cout << fixed;
	cout << scientific;

	ofstream out("outfile.dat");
	auto coutbuf = cout.rdbuf(out.rdbuf());

	string infile;
	string method;
	string basis_set;
	int sps;
	double eps;
	int max_cycles;
	string pop;

	cin >> infile;
	cin >> method;
	cin >> basis_set;
	cin >> sps;
	cin >> eps;
	cin >> max_cycles;
	cin >> pop;

	int N;
	int K;
	int nupdown;
	int Na;
	int Nb;
	bool r;
	double S2_e;
	double S2_UHF;
	int cycles = 1;
	double err = eps + 1;
	vector<vector<vector<vector<double>>>> eris;

	double nuc;
	double Eo;
	double temp_Eo;
	double E_tot;
	int icd = 0;

	//////////////////////////////////////////////////
	// s		// overlap			//
	// t		// kinetic			//
	// v		// nuclear attraction		//
	// hcore	// core Hamiltonian		//
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
	cout << "JESSE DICENSO\n";
	cout << "SUMMER 2024 - (probably May 2030?)\n\n";

	cout << "Reading input...\n\n";
	
	Molecule M(infile, basis_set);
	N = M.Nelec;
	K = M.AOs.size();
	nupdown = M.NUPDOWN;
	if(method=="RHF"){
		r = true;
	}
	else if(method=="UHF"){
		r = false;
	}
	if(r && (nupdown!=0)){
		cerr << "ERR: Restricted calculation cannot be performed for NUPDOWN!=0! (Try unrestricted...)\n";
		return 0;
	}
	Na = (N+nupdown)/2;
	Nb = (N-nupdown)/2;
	S2_e = 0.5*(Na-Nb)*(0.5*(Na-Nb)+1);
	cout << "\tinfile\t\t" << infile << '\n';
	cout << "\tmethod\t\t" << method << '\n';
	cout << "\tbasis\t\t"  << basis_set << '\n';
	cout << "\tsps\t\t\t" << sps << ' ';
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
	cout << "\teps\t\t\t" << setprecision(2) << eps << '\n';
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
	cout << "Success! There are " << N << " electrons " << "(" << Na << " alpha and " << Nb << " beta) and " << K << " basis functions.\n"; 
	nuc = nucrepl(M.Zvals, M.xyz);
	cout << "Nuclear Repulsion Energy = " << nuc << " Ha\n\n";
    cout.flush();
	
	cout << "Computing ERIs..." << endl;
	eris = ERIs(M.AOs);
	cout << "ERIs computed.\n" << endl;

	if(r){
		cout << "Performing restricted Hartree Fock...\n";	
		cout << "Using core Hamiltonian as initial guess to Fock matrix.\n\n";
	}
	else{
		cout << "Performing unrestricted Hartree Fock...\n";	
		cout << "Using core Hamiltonian as initial guess to alpha/beta Fock matrices.\n\n";
	}
	cout << "                   *************\n";
	cout << "                   * BEGIN SCF *\n";
	cout << "                   *************\n\n";
	cout << "----------------------------------------------------\n";
	cout << setw(3) << "N" << setw(20) << "E" <<  setw(20) << "err" << setw(10) << "step\n";
	cout << "----------------------------------------------------\n";
	cout.flush();	

	Matrix s = overlap(M.AOs);
	Matrix hcore = kinetic(M.AOs) + nuclear(M.AOs, M.Zvals, M.xyz);
	Matrix x = m_inv_sqrt(s);
	
	// Restricted
	if(r){
		vector<Matrix> temp_e_c;
		
		Matrix p (K, K);
		Matrix f (K, K);
		Matrix fo(K, K);
		Matrix e (K, K);
		Matrix co(K, K);
		Matrix c (K, K);
		
		f = R_F(hcore, p, eris);
		Eo = R_E0(p, hcore, f);
		fo = transpose(x) * f * x;
		temp_e_c = diagonalize(fo);
		e = temp_e_c[0];
		co = temp_e_c[1];
		c = x * co;
		p = R_density_matrix(c, N);

		if(sps==0){
			while((abs(err) > eps) && (cycles <= max_cycles)){
				R_FPI(s, hcore, eris, x, &p, &f, &fo, &e, &co, &c, &Eo, &err, N);
				cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err << setw(10) << "fp\n";
				cout.flush();
				cycles+=1;
			}
		}
		else{
			vector<Matrix> SPf;
			vector<Matrix> SPe;
			while((abs(err) > eps) && (cycles <= max_cycles)){
				R_DIIS(s, hcore, eris, x, &p, &f, &fo, &e, &co, &c, &Eo, &err, N, cycles, SPf, SPe, sps, &icd); // prev: SPf/e.data()
				if(icd!=0){
					break;
				}
				cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err;
				if(cycles <= 3){
					cout << setw(10) << "fp\n";
				}
				else if(cycles < sps){
					cout << setw(9) << "diis(" << SPe.size() << ")\n";
				}
				else{
					cout << setw(9) << "diis(" << sps << ")\n";
				}
				cout.flush();
				cycles+=1;
			}
			if(icd!=0){
			int fpi_forced_three;
				while((abs(err) > eps) && (cycles <= max_cycles) || (fpi_forced_three < 3)){
					R_FPI(s, hcore, eris, x, &p, &f, &fo, &e, &co, &c, &Eo, &err, N);
					cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err << setw(10) << "fp\n";
					cout.flush();
					fpi_forced_three+=1;
					cycles+=1;
				}
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

			cout << "Total E     = " << Eo + nuc << " Ha\n\n";
//
			grid mol_grid(M);
			double trps = Tr(p*s);
			double I_density = integrate_density(mol_grid, M, p);
			cout << "Trace of PS = " << trps << '\n';
			cout << "Integral of density = " << I_density << "\n\n";
//		
			R_print_orbitals(e, c, N, K);

			vector<double> popa(M.Zvals.size());
			if(pop=="lowdin"){
				popa = Lowdin_PA(M, p, s);
				cout << "=======================\n";
				cout << "Löwdin Pop. Analysis\n";
			}
			else if(pop=="mulliken"){
				popa = Mulliken_PA(M, p, s);
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
				cout << setw(3) << i+1 << setw(20) << popa[i] << '\n';
				sum_chg += popa[i];
			}
			cout << "=======================\n";
			cout << "Sum of atomic charges = " << sum_chg << "\n\n";
		}
	}
	// Unrestricted
	else{
		vector<Matrix> temp_e_c_a;
		vector<Matrix> temp_e_c_b;
		
		Matrix pt (K, K);
		Matrix pa (K, K);
		Matrix pb (K, K);
		Matrix fa (K, K);
		Matrix fb (K, K);
		Matrix fao(K, K);
		Matrix fbo(K, K);
		Matrix ea (K, K);
		Matrix eb (K, K);
		Matrix cao(K, K);
		Matrix cbo(K, K);
		Matrix ca (K, K);
		Matrix cb (K, K);
		
		fa = UR_F(hcore, pt, pa, eris);
		fb = UR_F(hcore, pt, pb, eris);
		Eo = UR_E0(pt, pa, pb, hcore, fa, fb);
		fao = transpose(x) * fa * x;
		fbo = transpose(x) * fb * x;
		temp_e_c_a = diagonalize(fao);
		temp_e_c_b = diagonalize(fbo);
		ea = temp_e_c_a[0];
		cao = temp_e_c_a[1];
		eb = temp_e_c_b[0];
		cbo = temp_e_c_b[1];
		ca = x * cao;
		cb = x * cbo;
		pa = UR_density_matrix(ca, Na);
		pb = UR_density_matrix(cb, Nb);
		pt = pa + pb;

		if(sps==0){
			while((abs(err) > eps) && (cycles <= max_cycles)){
				UR_FPI(s, hcore, eris, x, &pt, &pa, &pb, &fa, &fb, &fao, &fbo, &ea, &eb, &cao, &cbo, &ca, &cb, &Eo, &err, Na, Nb);
				cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err << setw(10) << "fp\n";
				cout.flush();
				cycles+=1;
			}
		}
		else{
			vector<Matrix> SPfa;
			vector<Matrix> SPfb;
			vector<Matrix> SPea;
			vector<Matrix> SPeb;
			while((abs(err) > eps) && (cycles <= max_cycles)){
				UR_DIIS(s, hcore, eris, x, &pt, &pa, &pb, &fa, &fb, &fao, &fbo, &ea, &eb, &cao, &cbo, &ca, &cb, &Eo, &err, Na, Nb, cycles, SPfa, SPfb, SPea, SPeb, sps, &icd);
				if(icd!=0){
					break;
				}
				cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err;
				if(cycles <= 3){
					cout << setw(10) << "fp\n";
				}
				else if(cycles < sps){
					cout << setw(9) << "diis(" << SPea.size() << ")\n";
				}
				else{
					cout << setw(9) << "diis(" << sps << ")\n";
				}
				cout.flush();
				cycles+=1;
			}
			if(icd!=0){
			int fpi_forced_three;
				while((abs(err) > eps) && (cycles <= max_cycles) || (fpi_forced_three < 3)){
					UR_FPI(s, hcore, eris, x, &pt, &pa, &pb, &fa, &fb, &fao, &fbo, &ea, &eb, &cao, &cbo, &ca, &cb, &Eo, &err, Na, Nb);
					cout << setw(3) << cycles << setw(20) << Eo << setw(20) << err << setw(10) << "fp\n";
					cout.flush();
					fpi_forced_three+=1;
					cycles+=1;
				}
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

			S2_UHF = S2_e + Nb;
			for(int i = 0; i < Na; i++){
				for(int j = 0; j < Nb; j++){
					double tempsum = 0;
					for(int mu = 0; mu < K; mu++){
						for(int nu = 0; nu < K; nu++){
							tempsum += (ca.matrix[mu][i])*(cb.matrix[nu][j])*(s.matrix[mu][nu]);
						}
					}
					S2_UHF -= tempsum*tempsum;
				}
			}

			cout << "Total E     = " << Eo + nuc << " Ha\n";
			cout << "Exact <S^2> = " << S2_e << '\n';
			cout << "UHF   <S^2> = " << S2_UHF << "\n\n";		
//
			grid mol_grid(M);
			double trpas = Tr(pa*s);
			double trpbs = Tr(pb*s);
			double trpts = Tr(pt*s);
			double I_density_a = integrate_density(mol_grid, M, pa);
			double I_density_b = integrate_density(mol_grid, M, pb);
			double I_density_t = integrate_density(mol_grid, M, pt);
			cout << "Trace    of Pa*S  = " << trpas << '\n';
			cout << "Integral of rho_a = " << I_density_a << '\n';
			cout << "Trace    of Pb*S  = " << trpbs << '\n';
			cout << "Integral of rho_b = " << I_density_b << '\n';
			cout << "Trace    of Pt*S  = " << trpts << '\n';
			cout << "Integral of rho_t = " << I_density_t << "\n\n";
//
			UR_print_orbitals(ea, eb, ca, cb, Na, Nb, K);

			vector<double> popa(M.Zvals.size());
			if(pop=="lowdin"){
				popa = Lowdin_PA(M, pt, s);
				cout << "=======================\n";
				cout << "Löwdin Pop. Analysis\n";
			}
			else if(pop=="mulliken"){
				popa = Mulliken_PA(M, pt, s);
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
				cout << setw(3) << i+1 << setw(20) << popa[i] << '\n';
				sum_chg += popa[i];
			}
			cout << "=======================\n";
			cout << "Sum of atomic charges = " << sum_chg << "\n\n";
		}
	}

	cout << "\n***************************************\n";
	cout <<   "*                                     *\n";
	cout <<   "*     Thank you, have a nice day!     *\n";
	cout <<   "*                                     *\n";
	cout <<   "***************************************\n\n";
}

vector<double> Lowdin_PA(Molecule M, Matrix P, Matrix S){
	vector<double> LPA(M.Zvals.size());
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

vector<double> Mulliken_PA(Molecule M, Matrix P, Matrix S){
	vector<double> MPA(M.Zvals.size());
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
	cout << "***************\n";
	cout << "* MO Energies *\n";
	cout << "***************\n\n";
	cout << "Occupied:\n";
	for(int i = 0; i < Nocc/2; i++){
		cout << setw(4) << 'E' << i+1 << setw(20) << E.matrix[i][i] << '\n';
	}
	if(Kb>(Nocc/2)){
		cout << "Virtual:\n";
		for(int i = Nocc/2; i < Kb; i++){
			cout << setw(4) << 'E' << i+1 << setw(20) << E.matrix[i][i] << '\n'; 
		}
	}
	cout << '\n';
	cout << "*******************\n";
	cout << "* MO Coefficients *\n";
	cout << "*******************\n\n";
	cout << "Occupied:\n";
	int count = 0;
	for(int i = 0; i < Nocc/2; i++){
		cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
		for(int j = 0; j < Kb; j++){
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
	if(Kb>(Nocc/2)){
		cout << "Virtual:\n";
		for(int i = Nocc/2; i < Kb; i++){
			cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
			for(int j = 0; j < Kb; j++){
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

void UR_print_orbitals(Matrix Ea, Matrix Eb, Matrix Ca, Matrix Cb, int Nocca, int Noccb, int Kb){
	cout << "***************\n";
	cout << "* MO Energies *\n";
	cout << "***************\n\n";
	cout << "Occupied (alpha):\n";
	for(int i = 0; i < Nocca; i++){
		cout << setw(4) << 'E' << i+1 << setw(20) << Ea.matrix[i][i] << '\n';
	}
	if(Kb>Nocca){
		cout << "Virtual (alpha):\n";
		for(int i = Nocca; i < Kb; i++){
			cout << setw(4) << 'E' << i+1 << setw(20) << Ea.matrix[i][i] << '\n'; 
		}
	}
	cout << '\n';
	cout << "Occupied (beta):\n";
	for(int i = 0; i < Noccb; i++){
		cout << setw(4) << 'E' << i+1 << setw(20) << Eb.matrix[i][i] << '\n';
	}
	if(Kb>Noccb){
		cout << "Virtual (beta):\n";
		for(int i = Noccb; i < Kb; i++){
			cout << setw(4) << 'E' << i+1 << setw(20) << Eb.matrix[i][i] << '\n'; 
		}
	}
	cout << '\n';
	cout << "*******************\n";
	cout << "* MO Coefficients *\n";
	cout << "*******************\n\n";
	cout << "Occupied (alpha):\n";
	int count = 0;
	for(int i = 0; i < Nocca; i++){
		cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
		for(int j = 0; j < Kb; j++){
			cout << Ca.matrix[j][i] << setw(20);
			count++;
			if(count>=5){
				cout << '\n' << setw(25);
				count = 0;
			}
		}
		count = 0;
		cout << '\n';
	}
	if(Kb>Nocca){
		cout << "Virtual (alpha):\n";
		for(int i = Nocca; i < Kb; i++){
			cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
			for(int j = 0; j < Kb; j++){
				cout << Ca.matrix[j][i] << setw(20);
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
	cout << "Occupied (beta):\n";
	count = 0;
	for(int i = 0; i < Noccb; i++){
		cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
		for(int j = 0; j < Kb; j++){
			cout << Cb.matrix[j][i] << setw(20);
			count++;
			if(count>=5){
				cout << '\n' << setw(25);
				count = 0;
			}
		}
		count = 0;
		cout << '\n';
	}
	if(Kb>Noccb){
		cout << "Virtual (beta):\n";
		for(int i = Noccb; i < Kb; i++){
			cout << setw(5) << "MO" << i+1 << '\n' << setw(25); 
			for(int j = 0; j < Kb; j++){
				cout << Cb.matrix[j][i] << setw(20);
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
