#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <xmmintrin.h>
#include "../timing.h"

union ieee754_QNAN
{
	uint64_t i;
	double f;
	ieee754_QNAN() : i(0x7FFFFFFFFFFFFFFF) {}
};

const ieee754_QNAN absMask;
static const __m128d abs4Mask = _mm_set1_pd(absMask.f);

static const __m128d zero = _mm_set1_pd(0.0);
static const __m128d one = _mm_set1_pd(1.0);
static const __m128d half = _mm_set1_pd(0.5);
static const __m128d neghalf = _mm_set1_pd(-0.5);
static const __m128d k1 = _mm_set1_pd(1.330274429);
static const __m128d k2 = _mm_set1_pd(-1.821255978);
static const __m128d k3 = _mm_set1_pd(1.781477937);
static const __m128d k4 = _mm_set1_pd(-0.356563782);
static const __m128d k5 = _mm_set1_pd(0.31938153);

static const __m128d mul = _mm_set1_pd(0.2316419);
static const __m128d Log10 = _mm_set1_pd(2.302585);
static const __m128d invSqrt2Pi = _mm_set1_pd(0.39894228040);

static inline __m128d _mm_blendv_pd(__m128d a, __m128d b, __m128d c) {
        return _mm_or_pd(_mm_andnot_pd(c, a), _mm_and_pd(c, b));
}



double d_cnd[2];
double d_body[2];

__attribute__ ((always_inline))
__m128d cnd(__m128d X) {
	__m128d L = _mm_and_pd(abs4Mask, X);
	__m128d k = _mm_div_pd(one, _mm_add_pd(one, _mm_mul_pd(mul, L)));
	__m128d w = _mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(k1, k), k2), k), k3), k), k4), k), k5), k);
	_mm_store_pd(d_cnd, _mm_mul_pd(_mm_mul_pd(L, L), neghalf));
	d_cnd[0] = exp(d_cnd[0]); d_cnd[1] = exp(d_cnd[1]);
	w = _mm_mul_pd(_mm_mul_pd(w, invSqrt2Pi), _mm_load_pd(d_cnd));
	__m128d mask = _mm_cmpgt_pd(X, zero);
	return _mm_blendv_pd(w, _mm_sub_pd(one, w), mask); 
}


__attribute__ ((always_inline))
__m128d body(__m128d S, __m128d X, __m128d T, __m128d r, __m128d v) {
	__m128d tmp = _mm_mul_pd(v, _mm_sqrt_pd(T));
	_mm_store_pd(d_body, _mm_div_pd(S, X));
	d_body[0] = log(d_body[0]); d_body[1] = log(d_body[1]);
	__m128d d1 = _mm_div_pd(_mm_add_pd(_mm_div_pd(_mm_load_pd(d_body), Log10), _mm_mul_pd(_mm_add_pd(r, _mm_mul_pd(v, _mm_mul_pd(v, half))), T)), tmp);
	__m128d d2 = _mm_sub_pd(d1, tmp);
	_mm_store_pd(d_body, _mm_mul_pd(_mm_sub_pd(zero, r), T));
	d_body[0] = exp(d_body[0]); d_body[1] = exp(d_body[1]);
	__m128d result = _mm_sub_pd(_mm_mul_pd(S, cnd(d1)), _mm_mul_pd(X, _mm_mul_pd(_mm_load_pd(d_body), cnd(d2))));
	return result;
}

__attribute__ ((noinline))
double run(double* S, double* X, double* T, double* r, double* v) {
	__m128d sum = _mm_set1_pd(0);
	for(int j = 0; j < LENGTH; j+=2) {
		sum = _mm_add_pd(sum, body(_mm_load_pd(S+j), _mm_load_pd(X+j), _mm_load_pd(T+j), _mm_load_pd(r+j), _mm_load_pd(v+j)));
	}
	return ((double*)&sum)[0] + ((double*)&sum)[1];
}

int main(int argc, char** argv) {
	double* S = (double*)((uint64_t)malloc(sizeof(double)*(LENGTH+3)) >> 5 << 5);
	double* X = (double*)((uint64_t)malloc(sizeof(double)*(LENGTH+3)) >> 5 << 5);
	double* T = (double*)((uint64_t)malloc(sizeof(double)*(LENGTH+3)) >> 5 << 5);
	double* r = (double*)((uint64_t)malloc(sizeof(double)*(LENGTH+3)) >> 5 << 5);
	double* v = (double*)((uint64_t)malloc(sizeof(double)*(LENGTH+3)) >> 5 << 5);
	
	for(int i = 0; i < LENGTH; i++) {
		S[i] = 100;
		X[i] = 98;
		T[i] = 2;
		r[i] = 0.02;
		v[i] = 5;
	}

	start_timing();
	double sum = 0;
	for(int i = 0; i < ROUNDS; i++) {
		sum += run(S, X, T, r, v);
	}
	printf("%f \t(%f)\n", end_timing(), sum / (LENGTH * ROUNDS));

	return 0;
}


