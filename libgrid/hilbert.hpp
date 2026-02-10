#ifndef HILBERTHEADERDEF
#define HILBERTHEADERDEF

#include <cassert>
#include <cstdint>

void axes_to_transpose(uint32_t* X, const int b=21, const int n=3);
uint64_t transpose_to_hilbert(const uint32_t* const X, const int b=21, const int n=3);

#endif
