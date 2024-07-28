#include "1e.hpp"

double T(PGF A, PGF B){
	double Ix;
	double Iy;
	double Iz;
	double Bexp = B.getexp();
	double Bl = B.getl();
	double Bm = B.getm();
	double Bn = B.getn();
	double BN = B.getN();
	double Sab = S(A, B);
	PGF Bup2x(Bexp, B.xyz, Bl+2, Bm, Bn, BN);
	PGF Bup2y(Bexp, B.xyz, Bl, Bm+2, Bn, BN);
	PGF Bup2z(Bexp, B.xyz, Bl, Bm, Bn+2, BN);
	// Ix
	if(Bl < 2){
		Ix = Bexp*(2*Bl+1)*Sab - 2*(Bexp*Bexp)*S(A,Bup2x);
	}
	else{
		PGF Bdn2x(Bexp, B.xyz, (Bl-2), Bm, Bn, BN);
		Ix = Bexp*(2*Bl+1)*Sab - 2*(Bexp*Bexp)*S(A,Bup2x) - (Bl*(Bl-1)/2)*S(A,Bdn2x);
	}
	// Iy
	if(Bm < 2){
		Iy = Bexp*(2*Bm+1)*Sab - 2*(Bexp*Bexp)*S(A,Bup2y);
	}
	else{
		PGF Bdn2y(Bexp, B.xyz, Bl, (Bm-2), Bn, BN);
		Iy = Bexp*(2*Bm+1)*Sab - 2*(Bexp*Bexp)*S(A,Bup2y) - (Bm*(Bm-1)/2)*S(A,Bdn2y);
	}
	// Iz
	if(Bn < 2){
		Iz = Bexp*(2*Bn+1)*Sab - 2*(Bexp*Bexp)*S(A,Bup2z);
	}
	else{
		PGF Bdn2z(Bexp, B.xyz, Bl, Bm, (Bn-2), BN);	
		Iz = Bexp*(2*Bn+1)*Sab - 2*(Bexp*Bexp)*S(A,Bup2z) - (Bn*(Bn-1)/2)*S(A,Bdn2z);
	}
	// Normalization included in overlaps
	return (Ix + Iy + Iz);
}

double T(CGF A, CGF B){
        double sum = 0;
        std::vector<double> dA = A.getd();
        std::vector<double> dB = B.getd();
        std::vector<PGF> GA = A.getpgf();
        std::vector<PGF> GB = B.getpgf();
        for(int i = 0; i < A.getlen(); i++){
                for(int j = 0; j < B.getlen(); j++){
                        sum += dA[i] * (T(GA[i], GB[j]) / (GA[i].getN() * GB[j].getN())) * dB[j];
                }
        }
        return (A.getN() * B.getN() * sum);
}

double V(PGF A, PGF B){

	return 0;
}

double V(CGF A, CGF B){

	return 0;
}
