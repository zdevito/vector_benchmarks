
#include <math.h>
#include <stdio.h>
#include "../timing.h"
#define N 65536
#define ROUNDS 20

__attribute__ ((noinline)) 
double cnd(double X) {
	double L = fabs(X);
	double k = 1.0/(1.0+0.2316419*L);
	double k2 = k*k;
	double k3 = k2*k;
	double k4 = k2*k2;
	double k5 = k3*k2;
	double invSqrt2Pi = 0.39894228040;
	double w = 0.31938153*k + -0.356563782*k2 + 1.781477937*k3 + -1.821255978*k4 + 1.330274429*k5;
	w = w*invSqrt2Pi * exp(-L * L * 0.5);
	return X > 0 ? 1.0-w : w;
}

__attribute__ ((noinline))
double body(double S, double X, double T, double r, double v) {
	double d1 = (log(S/X)/log(10) + (r+v*v*0.5)*T) / (v * sqrt(T));
	double d2 = d1 - v * sqrt(T);
	double result = S * cnd(d1) - X * exp(-r * T) * cnd(d2);
	return result;
}

int main(int argc, char** argv) {

	double sum = 0;
	
	for(int i = 0; i < ROUNDS; i++) {
		for(int j = 0; j < N; j++) {
			sum += body(j /* don't lift this out of the loop */, 98, 2, 0.02, 5);
		}
	}
	start_timing();
	for(int i = 0; i < ROUNDS; i++) {
		for(int j = 0; j < N; j++) {
			sum += body(j /* don't lift this out of the loop */, 98, 2, 0.02, 5);
		}
	}

	printf("%f (time %f)\n", sum / (N * ROUNDS),end_timing());

	return 0;
}

