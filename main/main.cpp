#include "main.hpp"

using namespace std;

int main(int argc, char* argv[]){
	cout << fixed;
	cout << scientific;

	//ofstream out("outfile.dat");
	//auto coutbuf = cout.rdbuf(out.rdbuf());

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

	int cycles = 1;
	double err = eps + 1;

	double Eo;
	double temp_Eo;
	double E_tot;
	int icd = 0;

	//////////////////////////////////////////////////
	// s		// overlap							//
	// t		// kinetic							//
	// v		// nuclear attraction				//
	// hcore	// core Hamiltonian					//
	// x		// orthogonalization (symmetric)	//
	// p		// density							//
	// c		// coefficient						//
	// c		// orthogonalized coefficient		//
	// f		// fock								//
	// fo		// orthogonalized fock				//
	// e		// orbital energies					//
	//////////////////////////////////////////////////
	
	cout << "*************************\n";
	cout << "*                       *\n";
	cout << "*   QUANTUM CHEMISTRY   *\n";
	cout << "*                       *\n";
	cout << "*************************\n\n";
	cout << "JESSE DICENSO\n";
	cout << "SUMMER 2024 - (probably May 2030?)\n\n";

	cout << "Reading input...\n\n";
		
	const Molecule M(infile, basis_set);
	const int N = M.Nelec;
	const int K = M.AOs.size();
	const int nupdown = M.NUPDOWN;
	const int Na = (N+nupdown)/2;
	const int Nb = (N-nupdown)/2;
	const double S2_e = 0.5*(Na-Nb)*(0.5*(Na-Nb)+1);

	bool r;
	if(method.substr(0,2)=="R_"){
		r = true;
	}
	else if(method.substr(0,2)=="U_"){
		r = false;
	}
	else{
		cerr << "ERR: Must specify method as restricted (R_...) or unrestricted (U_...)!\n";
		return 0;
	}
	if(r && (nupdown!=0)){
		cerr << "ERR: Restricted calculation forbidden for NUPDOWN!=0! (Try unrestricted...)\n";
		return 0;
	}

	cout << "\tinfile      " << infile << '\n';
	cout << "\tmethod      " << method << '\n';
	cout << "\tbasis       "  << basis_set << '\n';
	cout << "\tsps         " << sps << ' ';
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
	cout << "\teps         " << setprecision(2) << eps << '\n';
	cout << "\tmax_cycles  " << max_cycles << '\n';
	cout << "\tpopulation  " << pop << "\n\n";
	cout << setprecision(10);
	cout << "--------------------------------------------------------------\n";
	cout << "Z" << setw(20) << "x (Bohr)" << setw(20) << "y (Bohr)" << setw(20) <<  "z (Bohr)" << '\n';
	cout << "--------------------------------------------------------------\n";
	for(int i = 0; i < M.Zvals.size(); i++){
		cout << M.Zvals[i] << setw(20) << M.xyz[i][0] << setw(20) << M.xyz[i][1] << setw(20) << M.xyz[i][2] << '\n';
	}
	cout << "--------------------------------------------------------------\n";
	cout << "Success! There are " << N << " electrons " << "(" << Na << " alpha and " << Nb << " beta) and " << K << " basis functions.\n"; 
	const double nuc = nucrepl(M.Zvals, M.xyz);
	cout << "Nuclear Repulsion Energy = " << nuc << " Ha\n\n";
    cout.flush();
	
	cout << "Computing ERIs..." << endl;
	const vector<vector<vector<vector<double>>>> eris = ERIs(M.AOs);
	cout << "ERIs computed.\n" << endl;

	const grid mol_grid(M);
	XC xc(method);
	xc.eris = &eris;
	xc.g    = &mol_grid;
	xc.mol  = &M;

	if(r){
		cout << "Performing restricted calculation...\n";	
		cout << "Using core Hamiltonian as initial guess to Fock matrix.\n\n";
	}
	else{
		cout << "Performing unrestricted calculation...\n";	
		cout << "Using core Hamiltonian as initial guess to alpha/beta Fock matrices.\n\n";
	}
	cout << "                   *************\n";
	cout << "                   * BEGIN SCF *\n";
	cout << "                   *************\n\n";
	cout << "----------------------------------------------------\n";
	cout << setw(3) << "N" << setw(20) << "E" <<  setw(20) << "err" << setw(10) << "step\n";
	cout << "----------------------------------------------------\n";
	cout.flush();	

	const Matrix s = overlap(M.AOs);
	const Matrix hcore = kinetic(M.AOs) + nuclear(M.AOs, M.Zvals, M.xyz);
	const Matrix x = m_inv_sqrt(s);
	
	// Restricted
	if(r){
		vector<Matrix> temp_e_c;
		
		Matrix p  (K, K);
		Matrix fxc(K, K);
		
		xc.P   = &p;
		xc.FXC = &fxc;

		Matrix f (K, K);
		Matrix fo(K, K);
		Matrix e (K, K);
		Matrix co(K, K);
		Matrix c (K, K);
		
		f = hcore;
		Eo = E0(xc, hcore, zero(K,K));
		fo = transpose(x) * f * x;
		temp_e_c = diagonalize(fo);
		e = temp_e_c[0];
		co = temp_e_c[1];
		c = x * co;
		p = R_density_matrix(c, N);

		if(sps==0){
			while((abs(err) > eps) && (cycles <= max_cycles)){
				R_FPI(s, hcore, x, &xc, &f, &fo, &e, &co, &c, &Eo, &err, N, cycles);
				cycles+=1;
			}
		}
		else{
			vector<Matrix> SPf;
			vector<Matrix> SPe;
			while((abs(err) > eps) && (cycles <= max_cycles)){
				R_DIIS(s, hcore, x, &xc, &f, &fo, &e, &co, &c, &Eo, &err, N, cycles, SPf, SPe, sps, &icd);
				if(icd!=0){
					break;
				}
				cycles+=1;
			}
			if(icd!=0){
			int fpi_forced_three;
				while((abs(err) > eps) && (cycles <= max_cycles) || (fpi_forced_three < 3)){
					R_FPI(s, hcore, x, &xc, &f, &fo, &e, &co, &c, &Eo, &err, N, cycles);
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
			
			double trps = Tr(p*s);
			double I_density = integrate_quad(mol_grid, density2, M, p);
			cout << "Trace of PS = " << trps << '\n';
			cout << "Integral of density = " << I_density << "\n\n";
			
			R_print_orbital_energies(e, N, K);
			//R_print_orbitals(c, N, K);

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
		
		Matrix pt  (K, K);
		Matrix pa  (K, K);
		Matrix pb  (K, K);
		Matrix fxca(K, K);
		Matrix fxcb(K, K);
		
		xc.P     = &pt;
		xc.P_A   = &pa;
		xc.P_B   = &pb;
		xc.FXC_A = &fxca;	
		xc.FXC_B = &fxcb;	

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
		
		fa = hcore;
		fb = hcore;
		Eo = E0(xc, hcore, zero(K,K));
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
				UR_FPI(s, hcore, x, &xc, &fa, &fb, &fao, &fbo, &ea, &eb, &cao, &cbo, &ca, &cb, &Eo, &err, Na, Nb, cycles);
				cycles+=1;
			}
		}
		else{
			vector<Matrix> SPfa;
			vector<Matrix> SPfb;
			vector<Matrix> SPea;
			vector<Matrix> SPeb;
			while((abs(err) > eps) && (cycles <= max_cycles)){
				UR_DIIS(s, hcore, x, &xc, &fa, &fb, &fao, &fbo, &ea, &eb, &cao, &cbo, &ca, 
						&cb, &Eo, &err, Na, Nb, cycles, SPfa, SPfb, SPea, SPeb, sps, &icd);
				if(icd!=0){
					break;
				}
				cycles+=1;
			}
			if(icd!=0){
			int fpi_forced_three;
				while((abs(err) > eps) && (cycles <= max_cycles) || (fpi_forced_three < 3)){
					UR_FPI(s, hcore, x, &xc, &fa, &fb, &fao, &fbo, &ea, &eb, &cao, &cbo, &ca, &cb, &Eo, &err, Na, Nb, cycles);
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

			double S2_UHF = S2_e + Nb;
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
			cout << "Calc. <S^2> = " << S2_UHF << "\n\n";		
			
			double trpas = Tr(pa*s);
			double trpbs = Tr(pb*s);
			double trpts = Tr(pt*s);
			double I_density_a = integrate_quad(mol_grid, density2, M, pa);
			double I_density_b = integrate_quad(mol_grid, density2, M, pb);
			double I_density_t = integrate_quad(mol_grid, density2, M, pt);
			cout << "Trace    of Pa*S  = " << trpas << '\n';
			cout << "Integral of rho_a = " << I_density_a << '\n';
			cout << "Trace    of Pb*S  = " << trpbs << '\n';
			cout << "Integral of rho_b = " << I_density_b << '\n';
			cout << "Trace    of Pt*S  = " << trpts << '\n';
			cout << "Integral of rho_t = " << I_density_t << "\n\n";
			
			UR_print_orbital_energies(ea, eb, Na, Nb, K);
			//UR_print_orbitals(ca, cb, Na, Nb, K);

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
