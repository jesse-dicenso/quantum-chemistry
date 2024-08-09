#include "miscmath.hpp"

int fact(int n){
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

int dfact(int n){
	assert((n >= 0) || (n%2!=0));
	if(n==0){
		return 1;
	}
	else if(n==1){
		return 1;
	}
	else if(n > 0){
		return n*dfact(n-2);
	}
	else{
		return dfact(n+2) / (n+2);
	}
}

int binomial(int n, int k){
	assert((n >= 0) && (k >= 0) && (n >= k));
	return fact(n) / (fact(k) * fact(n - k));
}

double dot(std::vector<double> a, std::vector<double> b){
	assert(a.size()==b.size());
	double sum = 0;
	for(int i = 0; i < a.size(); i++){
		sum += a[i]*b[i];
	}
	return sum;
}
/*
double E(int i, int j, int t, double a, double b, double QAB){
	assert((i>=0) && (j>=0));
	double p = a + b;
	if((t < 0) || (t > (i+j))){
		return 0;
	}
	else if((t==0) || ((i+j)==0)){
		return exp(-a*b*QAB*QAB/p);
	}
	else if(i>0){
		return E(i-1, j, t-1, a, b, QAB)/(2*p) - (b*QAB/p)*E(i-1, j, t, a, b, QAB) + (t+1)*E(i-1, j, t+1, a, b, QAB);
	}
	else{
		return E(i, j-1, t-1, a, b, QAB)/(2*p) + (a*QAB/p)*E(i, j-1, t, a, b, QAB) + (t+1)*E(i, j-1, t+1, a, b, QAB);
	}
}
*/
