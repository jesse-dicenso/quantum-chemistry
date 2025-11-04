#include "quadrature.hpp"
#include <iostream>
#include <iomanip>

double f(double x, double y, double z){
	return exp(-sqrt(x*x + y*y + z*z));
}

int main(){
	double result = lebedev_gauss_chebyshev(f, 1, 100);
	std::cout << std::setprecision(16) << result << '\n';
	std::cout << std::setprecision(16) << 8*M_PI << '\n';
	std::cout << std::setprecision(16) << result-8*M_PI << '\n';
	return 0;
}
