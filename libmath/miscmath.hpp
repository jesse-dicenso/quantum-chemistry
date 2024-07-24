#ifndef MISCMATHHEADERDEF
#define MISCMATHHEADERDEF

#include<cmath>
#include<cassert>

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

#endif
