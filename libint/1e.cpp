#include "1e.hpp"

double T(PGF A, PGF B){
	double Ix;
	double Iy;
	double Iz;
	double Aexp = A.getexp();
	double Bexp = B.getexp();
	
	double Al = A.getl();
	double Am = A.getm();
	double An = A.getn();
	double AN = A.getN();
	double Bl = B.getl();
	double Bm = B.getm();
	double Bn = B.getn();
	double BN = B.getN();

	PGF Aup1x(Aexp, A.xyz, Al+1, Am, An, AN);
	PGF Aup1y(Aexp, A.xyz, Al, Am+1, An, AN);
	PGF Aup1z(Aexp, A.xyz, Al, Am, An+1, AN);
	PGF Bup1x(Bexp, B.xyz, Bl+1, Bm, Bn, BN);
	PGF Bup1y(Bexp, B.xyz, Bl, Bm+1, Bn, BN);
	PGF Bup1z(Bexp, B.xyz, Bl, Bm, Bn+1, BN);

	// Ix
	PGF Adn1x(Aexp, A.xyz, ((Al-1)+abs(Al-1))/2, Am, An, AN);
	PGF Bdn1x(Bexp, B.xyz, ((Bl-1)+abs(Bl-1))/2, Bm, Bn, BN);
	Ix = 0.5*Al*Bl*S(Adn1x,Bdn1x) + 2*Aexp*Bexp*S(Aup1x,Bup1x) - Aexp*Bl*S(Aup1x,Bdn1x) - Bexp*Al*S(Adn1x,Bup1x);
	// Iy
	PGF Adn1y(Aexp, A.xyz, Al, ((Am-1)+abs(Am-1))/2, An, AN);
	PGF Bdn1y(Bexp, B.xyz, Bl, ((Bm-1)+abs(Bm-1))/2, Bn, BN);	
	Iy = 0.5*Am*Bm*S(Adn1y,Bdn1y) + 2*Aexp*Bexp*S(Aup1y,Bup1y) - Aexp*Bm*S(Aup1y,Bdn1y) - Bexp*Am*S(Adn1y,Bup1y);

	// Iz
	PGF Adn1z(Aexp, A.xyz, Al, Am, ((An-1)+abs(An-1))/2, AN);	
	PGF Bdn1z(Bexp, B.xyz, Bl, Bm, ((Bn-1)+abs(Bn-1))/2, BN);	
	Iz = 0.5*An*Bn*S(Adn1z,Bdn1z) + 2*Aexp*Bexp*S(Aup1z,Bup1z) - Aexp*Bn*S(Aup1z,Bdn1z) - Bexp*An*S(Adn1z,Bup1z);

	// Normalization included in overlaps
	return (Ix + Iy + Iz);
}

// Asymmetric version of T; doesn't work!
/*
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
*/

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
