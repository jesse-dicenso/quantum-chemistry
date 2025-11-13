#include "density.hpp"

struct R_context{
	const Molecule 	*molecule;
	const Matrix	*Pmatrix;
	int 			atom_idx;
};

double R_atomic_density(double x, double y, double z, const Molecule &mol, const Matrix &P, int atom_i){
	assert(mol.NUPDOWN==0);
	double density = 0;
	double eval_gf;
	for(int i = 0; i < P.rows; i++){
	//	if(mol.AOs[i].atom_index == atom_i){
			eval_gf = mol.AOs[i].evaluate(x, y, z);
			for(int j = i+1; j < P.cols; j++){
	//			if(mol.AOs[j].atom_index == atom_i){
					density += 2 * P.matrix[i][j] * eval_gf * mol.AOs[j].evaluate(x, y, z);
	//			}
			}
			density += P.matrix[i][i] * eval_gf * eval_gf;
	//	}
	}
	return density;
}

double R_atomic_density_wrapper(double x, double y, double z, void* ctx){
	R_context* r_ctx = static_cast<R_context*>(ctx);
	if(r_ctx->molecule->Natoms == 1){
		return R_atomic_density(x, y, z, *r_ctx->molecule, *r_ctx->Pmatrix, 0);
	}
	return becke_weight(x, y, z, r_ctx->atom_idx, 3, *r_ctx->molecule) * 
		   R_atomic_density(x, y, z, *r_ctx->molecule, *r_ctx->Pmatrix, r_ctx->atom_idx);
}

double integrate_R(const Molecule &mol, const Matrix &P, int n){
	double result = 0;
	R_context ctx = {&mol, &P, 0};
	if(mol.Natoms == 1){
		result = lebedev_gauss_chebyshev(R_atomic_density_wrapper, &ctx, mol.xyz[0][0], mol.xyz[0][1], mol.xyz[0][2], 
										 bragg_slater_radii[mol.Zvals[0]-1], n);
	}
	else{
		for(int i = 0; i < mol.Natoms; i++){
			ctx.atom_idx = i;
			result += lebedev_gauss_chebyshev(R_atomic_density_wrapper, &ctx, mol.xyz[i][0], mol.xyz[i][1], mol.xyz[i][2], 
											  bragg_slater_radii[mol.Zvals[i]-1], n);
		}
	}
	return result;
}

// Corresponds to w_l in Becke paper; common k = 3
double becke_weight(double x, double y, double z, int l, int k, const Molecule &mol){
	assert(l < mol.Natoms);
	double r_i, r_j, R_ij, mu_ij, P_i, P_l = 0;
	double sum_P = 0;

	if(mol.heteronuclear){
	double u_ij, a_ij = 0;
		for(int i = 0; i < mol.Natoms; i++){
			P_i = 1;
			r_i = sqrt(	(x-mol.xyz[i][0])*(x-mol.xyz[i][0]) + 
						(y-mol.xyz[i][1])*(y-mol.xyz[i][1]) + 
						(z-mol.xyz[i][2])*(z-mol.xyz[i][2]));
			for(int j = 0; j < mol.Natoms; j++){
				if(j != i){
					r_j  = sqrt((x-mol.xyz[j][0])*(x-mol.xyz[j][0]) + 
								(y-mol.xyz[j][1])*(y-mol.xyz[j][1]) + 
								(z-mol.xyz[j][2])*(z-mol.xyz[j][2]));
					R_ij = sqrt((mol.xyz[i][0]-mol.xyz[j][0])*(mol.xyz[i][0]-mol.xyz[j][0]) + 
								(mol.xyz[i][1]-mol.xyz[j][1])*(mol.xyz[i][1]-mol.xyz[j][1]) + 
								(mol.xyz[i][2]-mol.xyz[j][2])*(mol.xyz[i][2]-mol.xyz[j][2]));
					mu_ij = (r_i - r_j) / R_ij;
					u_ij  = ( (bragg_slater_radii[mol.Zvals[i]] / bragg_slater_radii[mol.Zvals[j]]) - 1 ) / 
							( (bragg_slater_radii[mol.Zvals[i]] / bragg_slater_radii[mol.Zvals[j]]) + 1 );
					a_ij  = std::max(-0.5, std::min(0.5, u_ij / (u_ij * u_ij - 1) ));
					mu_ij += a_ij * ( 1 - mu_ij * mu_ij );
					P_i *= 0.5 * (1 - becke_step(mu_ij, k));
				}
			}
			if(i == l) {P_l = P_i;}
			sum_P += P_i;
		}
	}
	else{
		for(int i = 0; i < mol.Natoms; i++){
			P_i = 1;
			r_i = sqrt(	(x-mol.xyz[i][0])*(x-mol.xyz[i][0]) + 
						(y-mol.xyz[i][1])*(y-mol.xyz[i][1]) + 
						(z-mol.xyz[i][2])*(z-mol.xyz[i][2]));
			for(int j = 0; j < mol.Natoms; j++){
				if(j != i){
					r_j  = sqrt((x-mol.xyz[j][0])*(x-mol.xyz[j][0]) + 
								(y-mol.xyz[j][1])*(y-mol.xyz[j][1]) + 
								(z-mol.xyz[j][2])*(z-mol.xyz[j][2]));
					R_ij = sqrt((mol.xyz[i][0]-mol.xyz[j][0])*(mol.xyz[i][0]-mol.xyz[j][0]) + 
								(mol.xyz[i][1]-mol.xyz[j][1])*(mol.xyz[i][1]-mol.xyz[j][1]) + 
								(mol.xyz[i][2]-mol.xyz[j][2])*(mol.xyz[i][2]-mol.xyz[j][2]));
					mu_ij = (r_i - r_j) / R_ij;
					P_i *= 0.5 * (1 - becke_step(mu_ij, k));
				}
			}
			if(i == l) {P_l = P_i;}
			sum_P += P_i;
		}
	}
	return P_l / sum_P;
}

// "Continuous step function" corresponds to f_k in Becke paper
double becke_step(double mu, int k){
	assert(k >= 1);
	mu = std::max(-1.0, std::min(mu, 1.0));
	if(k==1){
		return 1.5 * mu - 0.5 * mu * mu * mu;
	}
	else{
		double temp = becke_step(mu, k-1);
		return 1.5 * temp - 0.5 * temp * temp * temp;
	}
}

const double bragg_slater_radii[118] = {
	1.3228082881354475,	// H
	2.645616576270895,	// He
	2.7401028825662843,	// Li
	1.9842124322031713,	// Be
	1.6062672070216149,	// B
	1.3228082881354475,	// C
	1.2283219818400586,	// N
	1.1338356755446692,	// O
	0.9448630629538911,	// F
	2.834589188861673,	// Ne
	3.401507026634008,	// Na
	2.834589188861673,	// Mg
	2.362157657384728,	// Al
	2.0786987384985607,	// Si
	1.8897261259077822,	// P
	1.8897261259077822,	// S
	1.8897261259077822,	// Cl
	3.401507026634008,	// Ar
	4.157397476997121,	// K
	3.401507026634008,	// Ca
	3.023561801452452,	// Sc
	2.645616576270895,	// Ti
	2.551130269975506,	// V
	2.645616576270895,	// Cr
	2.645616576270895,	// Mn
	2.645616576270895,	// Fe
	2.551130269975506,	// Co
	2.551130269975506,	// Ni
	2.551130269975506,	// Cu
	2.551130269975506,	// Zn
	2.456643963680117,	// Ga
	2.362157657384728,	// Ge
	2.1731850447939496,	// As
	2.1731850447939496,	// Se
	2.1731850447939496,	// Br
	3.590479639224786,	// Kr
	4.440856395883288,	// Rb
	3.7794522518155644,	// Sr
	3.401507026634008,	// Y
	2.9290754951570626,	// Zr
	2.7401028825662843,	// Nb
	2.7401028825662843,	// Mo
	2.551130269975506,	// Tc
	2.456643963680117,	// Ru
	2.551130269975506,	// Rh
	2.645616576270895,	// Pd
	3.023561801452452,	// Ag
	2.9290754951570626,	// Cd
	2.9290754951570626,	// In
	2.7401028825662843,	// Sn
	2.7401028825662843,	// Sb
	2.645616576270895,	// Te
	2.645616576270895,	// I
	3.9684248644063427,	// Xe
	4.913287927360234,	// Cs
	4.062911170701732,	// Ba
	3.684965945520175,	// La
	3.4959933329293973,	// Ce
	3.4959933329293973,	// Pr
	3.4959933329293973,	// Nd
	3.4959933329293973,	// Pm
	3.4959933329293973,	// Sm
	3.4959933329293973,	// Eu
	3.401507026634008,	// Gd
	3.307020720338619,	// Tb
	3.307020720338619,	// Dy
	3.307020720338619,	// Ho
	3.307020720338619,	// Er
	3.307020720338619,	// Tm
	3.307020720338619,	// Yb
	3.307020720338619,	// Lu
	2.9290754951570626,	// Hf
	2.7401028825662843,	// Ta
	2.551130269975506,	// W
	2.551130269975506,	// Re
	2.456643963680117,	// Os
	2.551130269975506,	// Ir
	2.551130269975506,	// Pt
	2.551130269975506,	// Au
	2.834589188861673,	// Hg
	3.590479639224786,	// Tl
	3.401507026634008,	// Pb
	3.023561801452452,	// Bi
	3.590479639224786,	// Po
	2.7401028825662843,	// At
	3.9684248644063427,	// Rn
	3.401507026634008,	// Fr
	4.062911170701732,	// Ra
	3.684965945520175,	// Ac
	3.401507026634008,	// Th
	3.401507026634008,	// Pa
	3.307020720338619,	// U
	3.307020720338619,	// Np
	3.307020720338619,	// Pu
	3.307020720338619,	// Am
	3.307020720338619,	// Cm
	3.307020720338619,	// Bk
	3.307020720338619,	// Cf
	3.307020720338619,	// Es
	3.307020720338619,	// Fm
	3.307020720338619,	// Md
	3.307020720338619,	// No
	3.307020720338619,	// Lr
	3.307020720338619,	// Rf
	3.307020720338619,	// Db
	3.307020720338619,	// Sg
	3.307020720338619,	// Bh
	3.307020720338619,	// Hs
	3.307020720338619,	// Mt
	3.307020720338619,	// Ds
	3.307020720338619,	// Rg
	3.307020720338619,	// Cn
	3.307020720338619,	// Nh
	3.307020720338619,	// Fl
	3.307020720338619,	// Mc
	3.307020720338619,	// Lv
	3.307020720338619,	// Ts
	3.307020720338619	// Og
};
