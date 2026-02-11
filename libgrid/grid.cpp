#include "grid.hpp"

grid::grid(const Molecule &mol, int num_radial, int num_angular, int becke_k){
	// radial  : Gauss-Chebyshev of the second kind
	// angular : lebedev (degree 230)
	assert((num_radial >= 0) && (num_angular == 230));
	num_gridpoints = num_radial * num_angular * mol.Natoms;
	x.resize(num_gridpoints);
	y.resize(num_gridpoints);
	z.resize(num_gridpoints);
	w.resize(num_gridpoints);
	double gauss_chebyshev_x, gauss_chebyshev_w, gauss_chebyshev_factor, r;
	int idx = 0;
    for(int atom = 0; atom < mol.Natoms; atom++){
        for(int i = 1; i <= num_radial; i++){
            gauss_chebyshev_x = std::cos(M_PI * i / (num_radial + 1));
            r = (bragg_slater_radii[mol.Zvals[atom] - 1] / 2) * (1 + gauss_chebyshev_x) / (1 - gauss_chebyshev_x);
            gauss_chebyshev_factor = r / sqrt(1 - gauss_chebyshev_x * gauss_chebyshev_x);
            gauss_chebyshev_factor *= gauss_chebyshev_factor * gauss_chebyshev_factor * 2;
            gauss_chebyshev_w = std::sin(M_PI * i / (num_radial + 1));
            gauss_chebyshev_w *= gauss_chebyshev_w * M_PI / (num_radial + 1) * gauss_chebyshev_factor;
            for(int j = 0; j < num_angular; j++){
                x[idx] = lebedev_x[j] * r + mol.xyz[atom][0];
                y[idx] = lebedev_y[j] * r + mol.xyz[atom][1];
                z[idx] = lebedev_z[j] * r + mol.xyz[atom][2];
                w[idx] = gauss_chebyshev_w * lebedev_w[j];
                if(mol.Natoms > 1){w[idx] *= becke_weight(mol, x[idx], y[idx], z[idx], atom, becke_k);}
                idx += 1;
            }
        }
    }
	assert(num_gridpoints == idx);
    hilbert_sort(x, y, z, w, num_gridpoints);
}

// Corresponds to w_me in Becke paper; common k = 3
double becke_weight(const Molecule &mol, double x, double y, double z, int atom_me, int k){
	assert(atom_me < mol.Natoms);
	double r_i, r_j, R_ij, mu_ij, P_i, P_me = 0;
	double sum_P = 0;
    for(int i = 0; i < mol.Natoms; i++){
        P_i = 1;
        r_i = sqrt( (x-mol.xyz[i][0])*(x-mol.xyz[i][0]) + 
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
                if(mol.heteronuclear){
                    const double u_ij  = ( (bragg_slater_radii[mol.Zvals[i]-1] / bragg_slater_radii[mol.Zvals[j]-1]) - 1 ) / 
                        ( (bragg_slater_radii[mol.Zvals[i]-1] / bragg_slater_radii[mol.Zvals[j]-1]) + 1 );
                    const double a_ij  = std::max(-0.5, std::min(0.5, u_ij / (u_ij * u_ij - 1) ));
                    mu_ij += a_ij * ( 1 - mu_ij * mu_ij );
                }
                P_i *= 0.5 * (1 - becke_step(mu_ij, k));
            }
        }
        if(i == atom_me) {P_me = P_i;}
        sum_P += P_i;
    }
	return P_me / sum_P;
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

// Sort gridpoints by their index on a Hilbert curve (default Hilbert depth of 21 in 3D to fit in a 64 bit unsigned int)
void hilbert_sort(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& w, const int size_g){
    // Transform coordinates: molecule -> [0,1]^3 -> (uint)[0,2^b)^3 -> Hilbert curve
    const double x_max = *std::max_element(x.begin(), x.end());
    const double y_max = *std::max_element(y.begin(), y.end());
    const double z_max = *std::max_element(z.begin(), z.end());

    const double x_min = *std::min_element(x.begin(), x.end());
    const double y_min = *std::min_element(y.begin(), y.end());
    const double z_min = *std::min_element(z.begin(), z.end());

    const double d_x = x_max - x_min;
    const double d_y = y_max - y_min;
    const double d_z = z_max - z_min;

    const double delta = ((d_x > d_y) ? ((d_x > d_z) ? d_x : d_z) : ((d_y > d_z) ? d_y : d_z));

    const uint32_t M = (1u << 21) - 1;
    uint32_t X[3] = {0, 0, 0};
    std::vector<uint64_t> H(size_g);
    for(int i = 0; i < size_g; i++){
        X[0] = (uint32_t)(M * ((x[i] - x_min) / delta));
        X[1] = (uint32_t)(M * ((y[i] - y_min) / delta));
        X[2] = (uint32_t)(M * ((z[i] - z_min) / delta));
        axes_to_transpose(X);
        H[i] = transpose_to_hilbert(X);
    }
    // sort x, y, z, w according to H
    sort_by_indices(H, x);
    sort_by_indices(H, y);
    sort_by_indices(H, z);
    sort_by_indices(H, w);
}

void sort_by_indices(const std::vector<uint64_t>& v1, std::vector<double>& v2){
    assert(v1.size() == v2.size());
    size_t size_v = v1.size();
    std::vector<size_t> idx(size_v);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&v1](size_t i1, size_t i2){return v1[i1] < v1[i2];});

    std::vector<double> v2_tmp(size_v);
    for(size_t i = 0; i < size_v; i++){
        v2_tmp[i] = v2[idx[i]];
    }
    v2 = v2_tmp;
}
