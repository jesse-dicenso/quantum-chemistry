#ifndef XCHEADERDEF
#define XCHEADERDEF

#include "dft_helper.hpp"

#include <string>
#include <functional>
#include <unordered_map>


class XC_inp{
	public:
		XC_inp(const std::string& method_name);
		
		std::string method;
	
		bool is_HF;
		bool is_LDA;
		bool is_GGA;
		
		Matrix* PT = nullptr;
		Matrix* PA = nullptr;
		Matrix* PB = nullptr;

		std::vector<std::vector<std::vector<std::vector<double>>>>* eris = nullptr;
		Molecule* mol = nullptr;
		grid* g = nullptr;
};

struct XC_ret{
	Matrix F_XC_1;	// F_XC if restricted; F_XC_a if unrestricted
	Matrix F_XC_2;	// F_XC_b if unrestricted
};

extern std::unordered_map<std::string, std::function<XC_ret(const XC_inp&)>> xc_register;

XC_ret F_XC(XC_inp* inp);

// HF //
XC_ret R_HF_X(const XC_inp& inp);
XC_ret U_HF_X(const XC_inp& inp);
// LDA //
XC_ret R_Slater_X(const XC_inp& inp);
XC_ret U_Slater_X(const XC_inp& inp);

#endif
