#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../timing.h"

__attribute__ ((always_inline))
double cnd(double X) {
	double L = fabs(X);
	double k = 1.0/(1.0+0.2316419*L);
	double w = (((((((1.330274429*k) - 1.821255978)*k) + 1.781477937)*k) - 0.356563782)*k + 0.31938153)*k;
	double invSqrt2Pi = 0.39894228040;
	w = w*invSqrt2Pi * exp(X * X * -0.5);
	return X > 0 ? 1.0-w : w;
}

double global = 0;

__attribute__ ((always_inline))
double body(double S, double X, double T, double r, double v) {
	double tmp = v * sqrt(T);
	double d1 = (log(S/X)/2.302585 + (r+v*v*0.5)*T) / tmp;
	double d2 = d1 - tmp;
	double result = S * cnd(d1) - X * exp(-r * T) * cnd(d2);
	global += result; // avoid some loop invariant code motion
	return result;
}

__attribute__ ((noinline))
double run(double* S, double* X, double* T, double* r, double* v) {
	double sum = 0;
	for(int j = 0; j < LENGTH; j++) {
		sum += body(S[j], X[j], T[j], r[j], v[j]);
	}
	return sum;
}

int main(int argc, char** argv) {

	double* S = (double*)((uint64_t)malloc(sizeof(double)*LENGTH+3) >> 4 << 4);
	double* X = (double*)((uint64_t)malloc(sizeof(double)*LENGTH+3) >> 4 << 4);
	double* T = (double*)((uint64_t)malloc(sizeof(double)*LENGTH+3) >> 4 << 4);
	double* r = (double*)((uint64_t)malloc(sizeof(double)*LENGTH+3) >> 4 << 4);
	double* v = (double*)((uint64_t)malloc(sizeof(double)*LENGTH+3) >> 4 << 4);
	
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

