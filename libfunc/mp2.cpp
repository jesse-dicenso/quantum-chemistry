#include "mp2.hpp"

struct tens4 {
    size_t m_dims[4];
    size_t m_strides[3];
    std::vector<double> m_data;

    tens4(size_t d1, size_t d2, size_t d3, size_t d4) : m_dims{d1, d2, d3, d4}, m_strides{d2*d3*d4, d3*d4, d4} {
        m_data.assign(d1*d2*d3*d4, 0.0);
    }
    const double& operator()(size_t i, size_t j, size_t k, size_t l) const {
        assert((i < m_dims[0]) && (j < m_dims[1]) && (k < m_dims[2]) && (l < m_dims[3]));
        return m_data[i*m_strides[0] + j*m_strides[1] + k*m_strides[2] + l];
    }
    double& operator()(size_t i, size_t j, size_t k, size_t l) {
        assert((i < m_dims[0]) && (j < m_dims[1]) && (k < m_dims[2]) && (l < m_dims[3]));
        return m_data[i*m_strides[0] + j*m_strides[1] + k*m_strides[2] + l];
    }
};

//
//  Restricted MP2 correlation energy
//
double MP2_ENERGY(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, // AO ERIs
        const Matrix& C, // MO coefficient matrix, K * K
        const Matrix& e, // Orbital energies (HF eigenvalues), K * K
        size_t N,        // # electrons
        size_t K         // # basis functions
    ) 
{  
    // 
    //  Perform standard O(N^5) AO->MO ERI transformation
    //
    const size_t nocc  = N/2;
    const size_t nvirt = K - nocc;

    // Transform first index (occ)
    tens4 T1(nocc, K, K, K);
    for (size_t i = 0; i < nocc; ++i) 
    for (size_t nu = 0; nu < K; ++nu) 
    for (size_t ld = 0; ld < K; ++ld)
    for (size_t sg = 0; sg < K; ++sg) {
        double t1 = 0.0;
        for (size_t mu = 0; mu < K; ++mu) {
            t1 += C(mu, i) * eris[mu][nu][ld][sg];
        }
        T1(i, nu, ld, sg) = t1;
    }

    // Transform second index (virt)
    tens4 T2(nocc, nvirt, K, K);
    for (size_t i = 0; i < nocc ; ++i)
    for (size_t a = 0; a < nvirt; ++a)
    for (size_t ld = 0; ld < K; ++ld)
    for (size_t sg = 0; sg < K; ++sg) {
        double t2 = 0.0;
        for (size_t nu = 0; nu < K; ++nu) {
            t2 += C(nu, a+nocc) * T1(i, nu, ld, sg); 
        }
        T2(i, a, ld, sg) = t2;
    }

    // Transform third index (occ)
    tens4 T3(nocc, nvirt, nocc, K);
    for (size_t i = 0; i < nocc ; ++i)
    for (size_t a = 0; a < nvirt; ++a)
    for (size_t j = 0; j < nocc ; ++j)
    for (size_t sg = 0; sg < K; ++sg) {
        double t3 = 0.0;
        for (size_t ld = 0; ld < K; ++ld) {
            t3 += C(ld, j) * T2(i, a, ld, sg); 
        }
        T3(i, a, j, sg) = t3;
    }

    // Transform fourth index (virt)
    tens4 MO_ERIS(nocc, nvirt, nocc, nvirt);
    for (size_t i = 0; i < nocc ; ++i)
    for (size_t a = 0; a < nvirt; ++a)
    for (size_t j = 0; j < nocc ; ++j)
    for (size_t b = 0; b < nvirt; ++b) {
        double mo_eris = 0.0;
        for (size_t sg = 0; sg < K; ++sg) {
            mo_eris += C(sg, b+nocc) * T3(i, a, j, sg);
        }
       MO_ERIS(i, a, j, b) = mo_eris; 
    }

    //
    //  Compute E_MP2
    //
    double E_MP2 = 0.0;
    for (size_t i = 0; i < nocc ; ++i)   // MO (occ)
    for (size_t a = 0; a < nvirt; ++a)   // MO (virt)
    for (size_t j = 0; j < nocc ; ++j)   // MO (occ)
    for (size_t b = 0; b < nvirt; ++b) { // MO (virt)
        const size_t A = a + nocc;
        const size_t B = b + nocc;
        const double Delta_ijab = e(i,i) + e(j,j) - e(A,A) - e(B,B);
        const double numerator  = MO_ERIS(i,a,j,b) * (2.0 * MO_ERIS(i,a,j,b) - MO_ERIS(i,b,j,a)); 
        E_MP2 += numerator / Delta_ijab;
    }

    return E_MP2;
}
