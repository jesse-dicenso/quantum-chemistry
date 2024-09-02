#include "miscmath.hpp"

double long fact(double long n){
	assert(n >= 0);
	if(n==0){
		return 1;
	}
	if(n==1){
		return 1;
	}
	else{
		return n*fact(n-1);
	}
}

double long dfact(double long n){
	assert((n >= 0) || ((int)n%2!=0));
	if(n==0){
		return 1;
	}
	else if(n==1){
		return 1;
	}
	else if(n>1){
		return n*dfact(n-2);
	}
	else{
		return dfact(n+2) / (n+2);
	}
}

double dot(std::vector<double> a, std::vector<double> b){
	assert(a.size()==b.size());
	double sum = 0;
	for(int i = 0; i < a.size(); i++){
		sum += a[i]*b[i];
	}
	return sum;
}

double boys(int n, double x){
	assert((x >= 0) && (n >= 0));
	if(x < 0.18){
		double sum = 0;
		for(int k = 0; k < 7; k++){
			sum += pow(-x, k)/(fact(k)*(2*n + 2*k + 1));
		}
		return sum;
	}
	else if(x > 19.35){
		return dfact(2*n-1)*sqrt(M_PI/pow(x, 2*n+1))/pow(2, n+1);
	}
	else{
		// J. Chem. Theory Comput. 2017, 13, 3636−3649
		double sum = 0;
		for(int k = 0; k < 51; k++){
			sum += dfact(2*n-1)*pow(2*x, k)/dfact(2*n+2*k+1);
		}
		return exp(-x)*sum;
	}
}
