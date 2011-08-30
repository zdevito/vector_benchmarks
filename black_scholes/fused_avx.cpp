
#include <math.h>
#include <stdio.h>
#include <xmmintrin.h>
#include "../timing.h"

// some defines first
union ieee754_QNAN
{
   const double f;
   struct
   {
      const uint64_t mantissa:52, exp:11, sign:1;
   };

   ieee754_QNAN() : f(0.0), mantissa(0xFFFFFFFFFFFFF), exp(0x7FF), sign(0x0) {}
};

const ieee754_QNAN absMask;
static const __m256d abs4Mask = _mm256_set1_pd(absMask.f);

static const __m256d zero = _mm256_set1_pd(0.0);
static const __m256d one = _mm256_set1_pd(1.0);
static const __m256d half = _mm256_set1_pd(-0.5);
static const __m256d k1 = _mm256_set1_pd(1.330274429);
static const __m256d k2 = _mm256_set1_pd(-1.821255978);
static const __m256d k3 = _mm256_set1_pd(1.781477937);
static const __m256d k4 = _mm256_set1_pd(-0.356563782);
static const __m256d k5 = _mm256_set1_pd(0.31938153);

static const __m256d mul = _mm256_set1_pd(0.2316419);
static const __m256d Log10 = _mm256_set1_pd(2.302585);
static const __m256d invSqrt2Pi = _mm256_set1_pd(0.39894228040);

__attribute__ ((always_inline))
__m256d cnd(__m256d X) {
	__m256d L = _mm256_and_pd(abs4Mask, X);
	__m256d k = __mm256_div_pd(one, _mm256_add_pd(one, _mm256_mul_pd(mul, L)));
	__m256d w = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(k1, k), k2), k), k3), k), k4), k), k5), k);
	w = _mm256_mul_pd(_mm256_mul_pd(L, L), half);
	((double*)&w)[0] = exp(((double*)&w)[0]);
	((double*)&w)[1] = exp(((double*)&w)[1]);
	w = _mm256_mul_pd(_mm256_mul_pd(w, invSqrt2Pi), w);
	__m256d mask = _mm256_cmpgt_pd(X, zero);
	return _mm256_blendv_pd(_mm256_sub_pd(one, w), w, mask); 
}

__attribute__ ((always_inline))
__m256d body(__m256d S, __m256d X, __m256d T, __m256d r, __m256d v) {
	__m256d tmp = _mm256_mul_pd(v, _mm256_sqrt_pd(T));
	__m256d w = _mm256_div_pd(S, X);
	((double*)&w)[0] = log(((double*)&w)[0]);
	((double*)&w)[1] = log(((double*)&w)[1]);
	__m256d d1 = _mm256_div_pd(_mm256_add_pd(_mm256_div_pd(w, Log10), _mm256_mul_pd(_mm256_add_pd(r, _mm256_mul_pd(v, _mm256_mul_pd(v, half))), T)), tmp);
	__m256d d2 = _mm256_sub_pd(d1, tmp);
	w = _mm256_mul_pd(_mm256_sub_pd(zero, r), T);
	((double*)&w)[0] = exp(((double*)&w)[0]);
	((double*)&w)[1] = exp(((double*)&w)[1]);
	__m256d result = _mm256_sub_pd(_mm256_mul_pd(S, cnd(d1)), _mm256_mul_pd(X, _mm256_mul_pd(w, cnd(d2))));
	return result;
}

__attribute__ ((noinline))
double run(double S, double X, double T, double r, double v) {
	__m256d sum = _mm256_set1_pd(0);
	__m256d _S = _mm256_set1_pd(S);
	__m256d _X = _mm256_set1_pd(X);
	__m256d _T = _mm256_set1_pd(T);
	__m256d _r = _mm256_set1_pd(r);
	__m256d _v = _mm256_set1_pd(v);
	for(int j = 0; j < LENGTH/4; j++) {
		// If we just use S, the compiler will lift this out of the loop.
		sum = _mm256_add_pd(sum, body(_mm256_set1_pd(j), _X, _T, _r, _v));
	}
	return ((double*)&sum)[0] + ((double*)&sum)[1];
}

int main(int argc, char** argv) {
	
	start_timing();
	double sum = 0;
	for(int i = 0; i < ROUNDS; i++) {
		sum += run(98, 2, 0.02, 5);
	}
	printf("%f \t(%f)\n", end_timing(), sum / (LENGTH * ROUNDS));

	return 0;
}

