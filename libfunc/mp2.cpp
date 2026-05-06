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
        assert((i < d1) && (j < d2) && (k < d3) && (l < d4));
        return m_data[i*m_strides[0] + j*m_strides[1] + k*m_strides[2] + l];
    }
}

//
//  Restricted MP2 correlation energy
//
double MP2_ENERGY(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& eris, // AO ERIs
        const Matrix& C, // MO coefficient matrix, K * K
        const Matrix& e, // Orbital energies (HF eigenvalues), K * K
        int N,           // # electrons
        int K            // # basis functions
    ) 
{  
    // 
    // Perform standard O(N^5) AO->MO ERI transformation
    //
    const size_t nocc  = N/2;
    const size_t nvirt = K - nocc;
    tens4 MO_ERIS(nocc, nvirt, nocc, nvirt);

    // Transform first index (occ)
    for (size_t i = 0; i < nocc; ++i) // MO (occ)
    for (size_t mu = 0; mu < K; ++mu) // AO
    for (size_t nu = 0; nu < K; ++nu) // AO
    for (size_t ld = 0; ld < K; ++ld) // AO
    for (size_t sg = 0; sg < K; ++sg) // AO
        MO_ERIS(i, nu, ld, sg) += C(mu, i) * eris[mu][nu][ld][sg]; 
    }

    // Transform second index (virt)
    for (size_t i = 0; i < nocc ; ++i) // MO (occ)
    for (size_t a = 0; a < nvirt; ++a) // MO (virt)
    for (size_t nu = 0; nu < K; ++nu)  // AO
    for (size_t ld = 0; ld < K; ++ld)  // AO
    for (size_t sg = 0; sg < K; ++sg)  // AO
        MO_ERIS(i, a, ld, sg) += C(nu, a+nocc) * MO_ERIS(i, nu, ld, sg); 
    }

    // Transform third index (occ)
    for (size_t i = 0; i < nocc ; ++i) // MO (occ)
    for (size_t a = 0; a < nvirt; ++a) // MO (virt)
    for (size_t j = 0; j < nocc ; ++j) // MO (occ)
    for (size_t ld = 0; ld < K; ++ld)  // AO
    for (size_t sg = 0; sg < K; ++sg)  // AO
        MO_ERIS(i, a, j, sg) += C(ld, j) * MO_ERIS(i, a, ld, sg); 
    }

    // Transform fourth index (virt)
    for (size_t i = 0; i < nocc ; ++i) // MO (occ)
    for (size_t a = 0; a < nvirt; ++a) // MO (virt)
    for (size_t j = 0; j < nocc ; ++j) // MO (occ)
    for (size_t b = 0; b < nvirt; ++b) // MO (virt)
    for (size_t sg = 0; sg < K; ++sg)  // AO
        MO_ERIS(i, a, j, b) += C(sg, b+nocc) * MO_ERIS(i, a, j, sg); 
    }

    // Compute E_MP2
    double E_MP2 = 0.0;
    for (size_t i = 0; i < nocc ; ++i)   // MO (occ)
    for (size_t a = 0; a < nvirt; ++a)   // MO (virt)
    for (size_t j = 0; j < nocc ; ++j)   // MO (occ)
    for (size_t b = 0; b < nvirt; ++b) { // MO (virt)
        const double Delta_ijab = e(i,i) + e(j,j) - e(a,a) - e(b,b);
        const double numerator  = MO_ERIS(i,a,j,b) * (2.0 * MO_ERIS(i,a,j,b) - MO_ERIS(i,b,j,a));
        E_MP2 += numerator / Delta_ijab;
    }

    return E_MP2;
}
