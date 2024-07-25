#ifndef GFHEADERDEF
#define GFHEADERDEF

#include <vector>

using namespace std;

class PGF{
	public:
		PGF(double exponent, vector<double> pos, int Ll, int Lm, int Ln);
		
		vector<double> xyz;

		double getexp();
		int getl();
		int getm();
		int getn();
		double getN();

		void printpgf();
		friend bool operator== (const PGF &p1, const PGF &p2);

	private:
		double exp;
		int l;
		int m;
		int n;
		double N;
};

class CGF{
	public:
		CGF(int lenC, vector<PGF> pgfC, vector<double> dC);
		
		int getlen();
		vector<PGF> getpgf();
		vector<double> getd();
		double getN();
		void printcgf();
		friend bool operator== (const CGF &c1, const CGF &c2);

	private:
		int len;
		vector<PGF> pgf;
		vector<double> d;
		double N;
};

// Functions for Gaussian Product Theorem

// fk: returns components of coefficients
// for contractions of arbitrary angular momentum
// in Gaussian Product Theorem
double fk(int k, int y1, int y2, double PA, double PB);

// Product Constant
double K(PGF p1, PGF p2);

// New Center
vector<double> P(PGF p1, PGF p2);

// overlap test
double overlap(PGF A, PGF B);
double overlap(CGF A, CGF B);

#endif
