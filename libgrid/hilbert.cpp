//#include "hilbert.hpp"
#include <cassert>
#include <cstdint>
#include <iostream>

// From Skilling 2004, Programming the Hilbert Curve
void axes_to_transpose(uint32_t* X, const int b, const int n) // position, # bits, dimension 
{
    const uint32_t M = 1 << (b - 1);
    uint32_t t;
    for(uint32_t Q = M; Q > 1; Q >>= 1){
        const uint32_t P = Q - 1;
        for(int i = 0; i < n; i++){
            if(X[i] & Q){X[0] ^= P;}
            else{
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
        }
    }
    for(int i = 0; i < n; i++){
        X[i] ^= X[i-1];
    }
    t = 0;
    for(uint32_t Q = M; Q > 1; Q >>= 1){
        if(X[n-1] & Q){t ^= Q - 1;}
    }
    for(int i = 0; i < n; i++){
        X[i] ^= t;
    }
}

uint64_t transpose_to_hilbert(const uint32_t* const X, const int b, const int n){
    assert(n * b <= 64);
    uint64_t H = 0;
    uint64_t t = 1;
    for(int j = 0; j < b; j++){
        for(int i = n-1; i >= 0; i--){
            if(X[i]>>j & 1){H |= t;}
            t <<= 1;
        }
    }
    return H;
}

int main(){
    uint32_t X[3] = {5, 10, 20};
    axes_to_transpose(X, 5, 3);
    uint64_t H = transpose_to_hilbert(X, 5, 3);
    std::cout << "H = " << H << " check 7865" << std::endl;
}
