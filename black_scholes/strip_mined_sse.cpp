#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <xmmintrin.h>
#include "../timing.h"

#define VW (BLOCK)

static inline __m128d _mm_blendv_pd(__m128d a, __m128d b, __m128d c) {
	return _mm_or_pd(_mm_andnot_pd(c, a), _mm_and_pd(c, b));
}

template<int N>
void neg_op(double const* i, double* o) {
	__m128d* xa = (__m128d*)i;
	__m128d xb = _mm_setzero_pd();
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_sub_pd(xb, xa[j]);
		//o[j] = -i[j];
	}
}

union ieee754_QNAN
{
	uint64_t i;
	double f;
	ieee754_QNAN() : i(0x7FFFFFFFFFFFFFFF) {}
};

const ieee754_QNAN absMask; 
static const __m128d abs4Mask = _mm_load1_pd( &absMask.f);

template<int N>
void abs_op(double const* i, double* o) {
	__m128d* xi = (__m128d*)i;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_and_pd( abs4Mask, xi[j]);
		//o[j] = fabs(i[j]);
	}
}

template<int N>
void log_op(double const* i, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = log(i[j]);
	}
}

template<int N>
void exp_op(double const* i, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = exp(i[j]);
	}
}

template<int N>
void sqrt_op(double const* i, double* o) {
	__m128d* xi = (__m128d*)i;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_sqrt_pd(xi[j]);
	}
}

template<int N>
void add_op(double const* a, double const* b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d* xb = (__m128d*)b;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_add_pd(xa[j], xb[j]);
		//o[j] = a[j] + b[j];
	}
}

template<int N>
void adds_op(double const* a, double b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d xb = _mm_set1_pd(b);
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_add_pd(xa[j], xb);
		//o[j] = a[j] + b;
	}
}

template<int N>
void sub_op(double const* a, double const* b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d* xb = (__m128d*)b;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_sub_pd(xa[j], xb[j]);
		//o[j] = a[j] - b[j];
	}
}

template<int N>
void mul_op(double const* a, double const* b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d* xb = (__m128d*)b;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_mul_pd(xa[j], xb[j]);
		//o[j] = a[j] * b[j];
	}
}

template<int N>
void muladd_op(double const* a, double const* b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d* xb = (__m128d*)b;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_add_pd(xo[j], _mm_mul_pd(xa[j], xb[j]));
		//o[j] = a[j] * b[j];
	}
}

template<int N>
void muls_op(double const* a, double b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d xb = _mm_set1_pd(b);
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_mul_pd(xa[j], xb);
		//o[j] = a[j] * b;
	}
}

template<int N>
void mulsadd_op(double const* a, double b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d xb = _mm_set1_pd(b);
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_add_pd(xo[j], _mm_mul_pd(xa[j], xb));
		//o[j] = a[j] * b;
	}
}

template<int N>
void addsmul_op(double a, double const* b, double* o) {
	__m128d xa = _mm_set1_pd(a);
	__m128d* xb = (__m128d*)b;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_mul_pd(xb[j], _mm_add_pd(xo[j], xa));
		//o[j] = a[j] * b;
	}
}

template<int N>
void muladds_op(double const* a, double b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d xb = _mm_set1_pd(b);
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_add_pd(xb, _mm_mul_pd(xa[j], xo[j]));
		//o[j] = a[j] * b;
	}
}

template<int N>
void div_op(double const* a, double const* b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d* xb = (__m128d*)b;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_div_pd(xa[j], xb[j]);
		//o[j] = a[j] / b[j];
	}
}

template<int N>
void divs_op(double const* a, double b, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d xb = _mm_set1_pd(b);
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_div_pd(xa[j], xb);
		//o[j] = a[j] / b[j];
	}
}

template<int N>
void rcp_op(double const* a, double* o) {
	__m128d* xa = (__m128d*)a;
	__m128d xb = _mm_set1_pd(1.0);
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_div_pd(xb, xa[j]);
		//o[j] = a[j] / b[j];
	}
}

template<int N>
void rep_op(double d, double* o) {
	__m128d* xo = (__m128d*)o;
	__m128d xb = _mm_set1_pd(d);
	for(int j = 0; j < N/2; j++) {
		xo[j] = xb;
		//o[j] = d;
	}
}

template<int N>
double sum_op(double const* i) {
	double s = 0;
	for(int j = 0; j < N; j++) {
		s += i[j];
	}
	return s;
}

template<int N>
void gt0_op(double const* i, double*o ) {
	__m128d* xa = (__m128d*)i;
	__m128d* xb = (__m128d*)o;
	__m128d zero = _mm_set1_pd(0.0);
	for(int j = 0; j < N/2; j++) {
		xb[j] = _mm_cmpgt_pd(xa[j], zero);
		//o[j] = i[j] > 0 ? 1.0 : 0.0;
	}
}

template<int N>
void sel_op(double const* s, double const* a, double const* b, double* o) {
	__m128d* xs = (__m128d*)s;
	__m128d* xa = (__m128d*)a;
	__m128d* xb = (__m128d*)b;
	__m128d* xo = (__m128d*)o;
	for(int j = 0; j < N/2; j++) {
		xo[j] = _mm_blendv_pd(xa[j], xb[j], xs[j]);
		//o[j] = s[j] > 0 ? a[j] : b[j];
	}
}

template<int N>
double* getV() {
	return new double[N];
}

double* t0 = getV<VW>();
double* t1 = getV<VW>();
double* t2 = getV<VW>();

double* cnd(double* X) {
	abs_op<VW>(X, t0);
	muls_op<VW>(t0, 0.2316419, t0);
	adds_op<VW>(t0, 1.0, t0);
	rcp_op<VW>(t0, t0);

	muls_op<VW>(t0, 1.330274429, t1);
	/*adds_op<VW>(t1, -1.821255978, t1);
	mul_op<VW>(t0, t1, t1);
	adds_op<VW>(t1, 1.781477937, t1);
	mul_op<VW>(t0, t1, t1);
	adds_op<VW>(t1, -0.356563782, t1);
	mul_op<VW>(t0, t1, t1);
	adds_op<VW>(t1, 0.31938153, t1);
	mul_op<VW>(t0, t1, t1);*/
	addsmul_op<VW>(-1.821255978, t0, t1);
	addsmul_op<VW>(1.781477937, t0, t1);
	addsmul_op<VW>(-0.356563782, t0, t1);
	addsmul_op<VW>(0.31938153, t0, t1);
	muls_op<VW>(t1, 0.39894228040, t1);

	mul_op<VW>(X, X, t0);
	muls_op<VW>(t0, -0.5, t0);
	exp_op<VW>(t0, t0);
	mul_op<VW>(t1, t0, t1);

	rep_op<VW>(1, t0);
	sub_op<VW>(t0, t1, t0);
	gt0_op<VW>(X, t2);
	sel_op<VW>(t2, t1, t0, t0);

	return t0;
}

double* s0 = getV<VW>();
double* s1 = getV<VW>();
double* s2 = getV<VW>();

double body(double const* S, double const* X, double const* T, double const* r, double const* v) {

	// log(S/X)/2.302585
	div_op<VW>(S, X, s0);
	log_op<VW>(s0, s0);
	divs_op<VW>(s0, 2.302585, s0);

	// (r+v*v*0.5)*T	
	mul_op<VW>(v, v, s1);
	muls_op<VW>(s1, 0.5, s1);
	add_op<VW>(s1, r, s1);

	// +
	muladd_op<VW>(s1, T, s0);

	// v * sqrt(T)	
	sqrt_op<VW>(T, s1);
	mul_op<VW>(v, s1, s1);

	// d1 = () / (v*sqrt(T))
	div_op<VW>(s0, s1, s0);

	// d2 = d1 - v*sqrt(T)
	sub_op<VW>(s0, s1, s1);

	// S * cnd(d1)
	mul_op<VW>(S, cnd(s0), s0);
	
	// exp(-r * T)
	neg_op<VW>(r, s2);
	mul_op<VW>(s2, T, s2);
	exp_op<VW>(s2, s2);

	// X * exp(-r * T) * cnd(d2)
	mul_op<VW>(X, s2, s2);
	mul_op<VW>(s2, cnd(s1), s2);

	sub_op<VW>(s0, s2, s0);	
	return sum_op<VW>(s0);
}

__attribute__ ((noinline))
double run(double* S, double* X, double* T, double* r, double* v) {
	double sum = 0;
	for(int j = 0; j < LENGTH; j+=BLOCK) {
		sum += body(S+j, X+j, T+j, r+j, v+j);
	}
	return sum;
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

