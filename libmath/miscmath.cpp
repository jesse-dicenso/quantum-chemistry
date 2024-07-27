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
	assert(n >= 0);
	if(n==0){
		return 1;
	}
	if(n==1){
		return 1;
	}
	else{
		return n*dfact(n-2);
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
