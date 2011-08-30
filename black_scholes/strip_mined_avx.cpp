#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <immintrin.h>
#include "../timing.h"

#define VW (BLOCK) 
#define VL ((LENGTH)/(BLOCK))

template<int N>
void neg_op(double* i, double* o) {
	__m256d* xa = (__m256d*)i;
	__m256d xb = _mm256_setzero_pd();
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_sub_pd(xb, xa[j]);
		//o[j] = -i[j];
	}
}

// some defines first
union ieee754_QNAN
{
   const double f;
   struct 
   {
      uint64_t mantissa:52, exp:11, sign:1; 
   };
   
   ieee754_QNAN() : f(0.0) {
	 mantissa = 0xFFFFFFFFFFFFF;
	exp = 0x7FF;
	sign = 0x0;
   } 
};

const ieee754_QNAN absMask; 
static const __m256d abs4Mask = _mm256_set1_pd( absMask.f);

template<int N>
void abs_op(double* i, double* o) {
	__m256d* xi = (__m256d*)i;
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_and_pd( abs4Mask, xi[j]);
		//o[j] = fabs(i[j]);
	}
}

template<int N>
void log_op(double* i, double* o) {
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_setzero_pd();
		//o[j] = log(i[j]);
	}
}

template<int N>
void exp_op(double* i, double* o) {
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_setzero_pd();
		//o[j] = exp(i[j]);
	}
}

template<int N>
void sqrt_op(double* i, double* o) {
	__m256d* xi = (__m256d*)i;
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_setzero_pd();
		//xo[j] = _mm256_sqrt_pd(xi[j]);
		//o[j] = sqrt(i[j]);
	}
}

template<int N>
void add_op(double* a, double* b, double* o) {
	__m256d* xa = (__m256d*)a;
	__m256d* xb = (__m256d*)b;
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_add_pd(xa[j], xb[j]);
		//o[j] = a[j] + b[j];
	}
}

template<int N>
void adds_op(double* a, double b, double* o) {
	__m256d* xa = (__m256d*)a;
	__m256d xb = _mm256_set1_pd(b);
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_add_pd(xa[j], xb);
		//o[j] = a[j] + b;
	}
}

template<int N>
void sub_op(double* a, double* b, double* o) {
	__m256d* xa = (__m256d*)a;
	__m256d* xb = (__m256d*)b;
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_sub_pd(xa[j], xb[j]);
		//o[j] = a[j] - b[j];
	}
}

template<int N>
void mul_op(double* a, double* b, double* o) {
	__m256d* xa = (__m256d*)a;
	__m256d* xb = (__m256d*)b;
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_mul_pd(xa[j], xb[j]);
		//o[j] = a[j] * b[j];
	}
}

template<int N>
void muls_op(double* a, double b, double* o) {
	__m256d* xa = (__m256d*)a;
	__m256d xb = _mm256_set1_pd(b);
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_mul_pd(xa[j], xb);
		//o[j] = a[j] * b;
	}
}

template<int N>
void div_op(double* a, double* b, double* o) {
	__m256d* xa = (__m256d*)a;
	__m256d* xb = (__m256d*)b;
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_setzero_pd();
		//xo[j] = _mm256_div_pd(xa[j], xb[j]);
		//o[j] = a[j] / b[j];
	}
}

template<int N>
void rep_op(double d, double* o) {
	__m256d* xo = (__m256d*)o;
	__m256d xb = _mm256_set1_pd(d);
	for(int j = 0; j < N/4; j++) {
		xo[j] = xb;
		//o[j] = d;
	}
}

template<int N>
double sum_op(double* i) {
	double s = 0;
	for(int j = 0; j < N; j++) {
		s += i[j];
	}
	return s;
}

template<int N>
void gt0_op(double* i, double*o ) {
	__m256d* xa = (__m256d*)i;
	__m256d* xb = (__m256d*)o;
	__m256d zero = _mm256_set1_pd(0.0);
	for(int j = 0; j < N/4; j++) {
		xb[j] = _mm256_cmpgt_pd(xa[j], zero);
		//o[j] = i[j] > 0 ? 1.0 : 0.0;
	}
}

template<int N>
void sel_op(double* s, double* a, double* b, double* o) {
	__m256d* xs = (__m256d*)s;
	__m256d* xa = (__m256d*)a;
	__m256d* xb = (__m256d*)b;
	__m256d* xo = (__m256d*)o;
	for(int j = 0; j < N/4; j++) {
		xo[j] = _mm256_blendv_pd(xa[j], xb[j], xs[j]);
		//o[j] = s[j] > 0 ? a[j] : b[j];
	}
}

template<int N>
double* getV() {
	double* d = new double[N+1];
	return (double*)(((int64_t)(d+1))>>5<<5);
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

double* S = getV<VW>();
double* X = getV<VW>();
double* T = getV<VW>();
double* r = getV<VW>();
double* v = getV<VW>();
double* s0 = getV<VW>();
double* s1 = getV<VW>();
double* s2 = getV<VW>();

double global = 0;
	
double body(double _S, double _X, double _T, double _r, double _v) {
	global++;			// Avoid LICM
	rep_op<VW>(_S, S);
	rep_op<VW>(_X, X);
	rep_op<VW>(_T, T);
	rep_op<VW>(_r, r);
	rep_op<VW>(_v, v);

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
double run(double S, double X, double T, double r, double v) {
	double sum = 0;
	for(int j = 0; j < LENGTH; j++) {
		sum += body(S, X, T, r, v);
	}
	return sum;
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

