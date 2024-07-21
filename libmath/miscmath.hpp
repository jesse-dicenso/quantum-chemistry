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

#endif
