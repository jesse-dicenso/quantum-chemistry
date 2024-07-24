#ifndef ZEROEHEADERDEF
#define ZEROEHEADERDEF

#include "../libgfs/gf.hpp"
#include "../libmath/miscmath.hpp"

double overlap(PGF A, PGF B){
	double S = 1.0;
	if(A==B){
		return S;
	}

	double l1 = A.getl();
	double l2 = B.getl();
	double m1 = A.getm();
	double m2 = B.getm();
	double n1 = A.getn();
	double n2 = B.getn();

	double Kab = K(A, B);
	double gam = A.getexp() + B.getexp();
	vector<double> Pab = P(A, B);
	vector<double> PA(3);
	vector<double> PB(3);

	PA[0] = Pab[0]-A.xyz[0];
	PA[1] = Pab[1]-A.xyz[1];
	PA[2] = Pab[2]-A.xyz[2];

	PB[0] = Pab[0]-B.xyz[0];
	PB[1] = Pab[1]-B.xyz[1];
	PB[2] = Pab[2]-B.xyz[2];
	
	double Ix = 0;
	double Iy = 0;
	double Iz = 0;
	
	for(int i = 0; i <= ((l1+l2))/2; i++){
		Ix += fk(2*i, l1, l2, PA[0], PB[0])*(dfact(2*i-1)/pow(2*gam,i))*pow((M_PI/gam),0.5);
	}	
	for(int i = 0; i <= ((m1+m2))/2; i++){
		Iy += fk(2*i, m1, m2, PA[1], PB[1])*(dfact(2*i-1)/pow(2*gam,i))*pow((M_PI/gam),0.5);
	}
	for(int i = 0; i <= ((n1+n2))/2; i++){
		Iz += fk(2*i, n1, n2, PA[2], PB[2])*(dfact(2*i-1)/pow(2*gam,i))*pow((M_PI/gam),0.5);
	}
	
	// Remember to normalize
	S = (A.getN() * B.getN())*Kab*Ix*Iy*Iz;
	assert(S <= 1.0);
	return S;
}
/*
double overlap(CGF A, CGF B){
	return 0;
}
*/

#endif
