#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include "../timing.h"

template<typename T>
T* malloc_aligned(uint64_t n, uint64_t log2_alignment) {
	void* a = malloc((sizeof(T))*(n+1));
	return (double*)((((uint64_t)a)+(sizeof(T)-1)) >> log2_alignment << log2_alignment);
}
