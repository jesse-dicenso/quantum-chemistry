#include "0e.hpp"

double overlap(PGF A, PGF B){
        if(A==B){
                return 1.0;
        }
        double l1 = A.getl();
        double l2 = B.getl();
        double m1 = A.getm();
        double m2 = B.getm();
        double n1 = A.getn();
        double n2 = B.getn();

        double Kab = K(A, B);
        double gam = A.getexp() + B.getexp();
        std::vector<double> Pab = P(A, B);
        std::vector<double> PA(3);
        std::vector<double> PB(3);

        PA[0] = Pab[0]-A.xyz[0];
        PA[1] = Pab[1]-A.xyz[1];
        PA[2] = Pab[2]-A.xyz[2];

        PB[0] = Pab[0]-B.xyz[0];
        PB[1] = Pab[1]-B.xyz[1];
        PB[2] = Pab[2]-B.xyz[2];

        double Ix = fk(0, l1, l2, PA[0], PB[0])*pow((M_PI/gam),0.5);
        double Iy = fk(0, m1, m2, PA[0], PB[0])*pow((M_PI/gam),0.5);
        double Iz = fk(0, n1, n2, PA[0], PB[0])*pow((M_PI/gam),0.5);

        for(int i = 1; i <= ((l1+l2))/2; i++){
                Ix += fk(2*i, l1, l2, PA[0], PB[0])*(dfact(2*i-1)/pow(2*gam,i))*pow((M_PI/gam),0.5);
        }
        for(int i = 1; i <= ((m1+m2))/2; i++){
                Iy += fk(2*i, m1, m2, PA[1], PB[1])*(dfact(2*i-1)/pow(2*gam,i))*pow((M_PI/gam),0.5);
        }
        for(int i = 1; i <= ((n1+n2))/2; i++){
                Iz += fk(2*i, n1, n2, PA[2], PB[2])*(dfact(2*i-1)/pow(2*gam,i))*pow((M_PI/gam),0.5);
        }
        // Remember to normalize
        return (A.getN() * B.getN())*Kab*Ix*Iy*Iz;
}

double overlap(CGF A, CGF B){
        if(A==B){
                return 1.0;
        }
        double sum = 0;
        std::vector<double> dA = A.getd();
        std::vector<double> dB = B.getd();
        std::vector<PGF> GA = A.getpgf();
        std::vector<PGF> GB = B.getpgf();
        for(int i = 0; i < A.getlen(); i++){
                for(int j = 0; j < B.getlen(); j++){
                        sum += dA[i] * (overlap(GA[i], GB[j]) / (GA[i].getN() * GB[j].getN())) * dB[j];
                }
        }
        //Remember to normalize; we "unnormalized" the PGFs above
        return (A.getN() * B.getN() * sum);
}
