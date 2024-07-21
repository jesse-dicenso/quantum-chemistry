#ifndef GFHEADERDEF
#define GFHEADERDEF

#include<vector>

using namespace std;

class PGF{
	public:
		PGF(double exponent, double Rx, double Ry, double Rz, int Ll, int Lm, int Ln);
		
		double x;
		double y;
		double z;

		int getexp();
		int getl();
		int getm();
		int getn();
		double getN();

		void printpgf();

	private:
		double exp;
		int l;
		int m;
		int n;
		double N;
};

class CGF{
	public:
		CGF(int len, vector<PGF> pgf, vector<double> d);
		
		double getN();
		void printcgf();

	private:
		int len;
		vector<PGF> pgf;
		vector<double> d;
		double N;
}

#endif
