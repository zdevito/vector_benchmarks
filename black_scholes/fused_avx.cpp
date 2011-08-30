#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <immintrin.h>
#include "../timing.h"

union ieee754_QNAN
{
	uint64_t i;
	double f;
	ieee754_QNAN() : i(0x7FFFFFFFFFFFFFFF) {}
};

const ieee754_QNAN absMask;
static const __m256d abs4Mask = _mm256_set1_pd(absMask.f);

static const __m256d zero = _mm256_set1_pd(0.0);
static const __m256d one = _mm256_set1_pd(1.0);
static const __m256d half = _mm256_set1_pd(0.5);
static const __m256d neghalf = _mm256_set1_pd(-0.5);
static const __m256d k1 = _mm256_set1_pd(1.330274429);
static const __m256d k2 = _mm256_set1_pd(-1.821255978);
static const __m256d k3 = _mm256_set1_pd(1.781477937);
static const __m256d k4 = _mm256_set1_pd(-0.356563782);
static const __m256d k5 = _mm256_set1_pd(0.31938153);

static const __m256d mul = _mm256_set1_pd(0.2316419);
static const __m256d Log10 = _mm256_set1_pd(2.302585);
static const __m256d invSqrt2Pi = _mm256_set1_pd(0.39894228040);

static inline __m256d _mm256_blendv_pd(__m256d a, __m256d b, __m256d c) {
        return _mm256_or_pd(_mm256_andnot_pd(c, a), _mm256_and_pd(c, b));
}

double * d_cnd;
double * d_body;
__attribute__ ((always_inline))
__m256d cnd(__m256d X) {
 
	__m256d L = _mm256_and_pd(abs4Mask, X);
	__m256d k = _mm256_div_pd(one, _mm256_add_pd(one, _mm256_mul_pd(mul, L)));
	__m256d w = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(k1, k), k2), k), k3), k), k4), k), k5), k);
	_mm256_store_pd(d_cnd, _mm256_mul_pd(_mm256_mul_pd(L, L), neghalf));
	d_cnd[0] = exp(d_cnd[0]); d_cnd[1] = exp(d_cnd[1]); d_cnd[2] = exp(d_cnd[2]); d_cnd[3] = exp(d_cnd[3]);
	w = _mm256_mul_pd(_mm256_mul_pd(w, invSqrt2Pi), _mm256_load_pd(d_cnd));
	__m256d mask = _mm256_cmpgt_pd(X, zero);
	return _mm256_blendv_pd(w, _mm256_sub_pd(one, w), mask); 
}

__attribute__ ((always_inline))
__m256d body(__m256d S, __m256d X, __m256d T, __m256d r, __m256d v) {
    
	__m256d tmp = _mm256_mul_pd(v, _mm256_sqrt_pd(T));
	_mm256_store_pd(d_body, _mm256_div_pd(S, X));
	d_body[0] = log(d_body[0]); d_body[1] = log(d_body[1]); d_body[2] = log(d_body[2]); d_body[3] = log(d_body[3]);
	__m256d d1 = _mm256_div_pd(_mm256_add_pd(_mm256_div_pd(_mm256_load_pd(d_body), Log10), _mm256_mul_pd(_mm256_add_pd(r, _mm256_mul_pd(v, _mm256_mul_pd(v, half))), T)), tmp);
	__m256d d2 = _mm256_sub_pd(d1, tmp);
	_mm256_store_pd(d_body, _mm256_mul_pd(_mm256_sub_pd(zero, r), T));
	d_body[0] = exp(d_body[0]); d_body[1] = exp(d_body[1]); d_body[2] = exp(d_body[2]); d_body[3] = exp(d_body[3]);
	__m256d result = _mm256_sub_pd(_mm256_mul_pd(S, cnd(d1)), _mm256_mul_pd(X, _mm256_mul_pd(_mm256_load_pd(d_body), cnd(d2))));
	return result;
}

__attribute__ ((noinline))
double run(double* S, double* X, double* T, double* r, double* v) {
	__m256d sum = _mm256_set1_pd(0);
	for(int j = 0; j < LENGTH; j+=4) {
		sum = _mm256_add_pd(sum, body(_mm256_load_pd(S+j), _mm256_load_pd(X+j), _mm256_load_pd(T+j), _mm256_load_pd(r+j), _mm256_load_pd(v+j)));
	}
	return ((double*)&sum)[0] + ((double*)&sum)[1] + ((double*)&sum)[2] + ((double*)&sum)[3];
}

int main(int argc, char** argv) {
  d_cnd = (double*)((uint64_t)malloc(sizeof(double)*12) >> 5 << 5);
  d_body = (double*) ((uint64_t) malloc(sizeof(double)*12) >> 5 << 5);
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

