
#include <iostream>
#include <stdint.h>
#include <math.h>

#define LONGVECTOR_LENGTH 256*256
#define VECTOR_WIDTH 256
#define ROUNDS 20

#define NOINLINE __attribute__ ((noinline))
//#define NOINLINE
#define APPLY(expr) for(int64_t i = 0; i < VECTOR_WIDTH; i++) { expr; }

NOINLINE void neg_d(double*__restrict__ a) { APPLY(a[i] = -a[i]) }
NOINLINE void abs_d(double*__restrict__ a) { APPLY(a[i] = fabs(a[i])) }
NOINLINE void log_d(double*__restrict__ a) { APPLY(a[i] = log(a[i])) }
NOINLINE void exp_d(double*__restrict__ a) { APPLY(a[i] = exp(a[i])) }
NOINLINE void sqrt_d(double*__restrict__ a) { APPLY(a[i] = sqrt(a[i])) }

NOINLINE void copy_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] = b[i]) }
NOINLINE void copy_d(double*__restrict__ a, double b) { APPLY(a[i] = b) }

NOINLINE void add_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] += b[i]) }
NOINLINE void adds_d(double*__restrict__ a, double b) { APPLY(a[i] += b) }

NOINLINE void sub_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] -= b[i]) }
NOINLINE void subs_d(double*__restrict__ a, double b) { APPLY(a[i] -= b) }
NOINLINE void rsub_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] = b[i]-a[i]) }
NOINLINE void rsubs_d(double*__restrict__ a, double b) { APPLY(a[i] = b-a[i]) }

NOINLINE void mul_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] = a[i]*b[i]) }
NOINLINE void muls_d(double*__restrict__ a, double b) { APPLY(a[i] *= b) }

NOINLINE void div_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] /= b[i]) }
NOINLINE void divs_d(double*__restrict__ a, double b) { APPLY(a[i] /= b) }

NOINLINE void rdiv_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] = b[i]/a[i]) }
NOINLINE void rdivs_d(double*__restrict__ a, double b) { APPLY(a[i] = b/a[i]) }

NOINLINE void sum_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[0] += b[i]) }
NOINLINE void sum_d(double*__restrict__ a, double*__restrict__ b, int64_t*__restrict__ group) { APPLY(a[group[i]] += b[i]) }

NOINLINE void selgt0_d(double*__restrict__ a, double*__restrict__ b, double*__restrict__ mask) 
	{ APPLY(a[i] = mask[i] > 0 ? a[i] : b[i]) }

double  d_reg[32][VECTOR_WIDTH] __attribute__((aligned (16)));
int64_t i_reg[32][VECTOR_WIDTH];

bool mask[VECTOR_WIDTH];

NOINLINE double* cnd(double* x) {
	double* const& d0 = d_reg[0];
	double* const& d1 = d_reg[1];
	
	copy_d(d0, x);
	abs_d(d0);
	muls_d(d0, 0.2316419);
	adds_d(d0, 1.0);
	copy_d(d1, 1.0);
	rdiv_d(d0, d1);

	copy_d(d1, 1.330274429);
	mul_d(d1, d0);
	adds_d(d1, -1.821255978);
	mul_d(d1, d0);
	adds_d(d1, 1.781477937);
	mul_d(d1, d0);
	adds_d(d1, -0.356563782);
	mul_d(d1, d0);
	adds_d(d1, 0.31938153);
	mul_d(d1, d0);

	muls_d(d1, 0.39894228040);
	copy_d(d0, x);
	neg_d(d0);
	mul_d(d0, x);
	muls_d(d0, 0.5);
	exp_d(d0);

	mul_d(d1, d0);

	copy_d(d0, d1);
	rsubs_d(d0, 1.0);	

	selgt0_d(d0, d1, x);
	return d0;
}


NOINLINE double body(double* S, double* X, double* T, double* r, double* v, double* result) {
	double* const& d0 = d_reg[2];
	double* const& d1 = d_reg[3];
	double* const& d2 = d_reg[4];
	
	copy_d(d0, S);
	div_d(d0, X);
	log_d(d0);
	divs_d(d0, 2.302585);
	
	copy_d(d1, v);
	mul_d(d1, v);
	muls_d(d1, 0.5);
	add_d(d1, r);
	mul_d(d1, T);
	
	add_d(d0, d1);

	
	copy_d(d1, T);
	sqrt_d(d1);
	mul_d(d1, v);
	
	div_d(d0, d1);

	rsub_d(d1, d0);

	copy_d(d0, cnd(d0));
	mul_d(d0, S);
	
	copy_d(d2, r);
	mul_d(d2, T);
	neg_d(d2);
	exp_d(d2);
	mul_d(d2, X);
	mul_d(d2, cnd(d1));

	sub_d(d0, d2);

	sum_d(result, d0);
}


int main(int argc, char** argv) {

	double* S = new double[LONGVECTOR_LENGTH];
	double* X = new double[LONGVECTOR_LENGTH];
	double* T = new double[LONGVECTOR_LENGTH];
	double* r = new double[LONGVECTOR_LENGTH];
	double* v = new double[LONGVECTOR_LENGTH];

	for(int64_t i = 0; i < LONGVECTOR_LENGTH; i++) {
		S[i] = 100;
		X[i] = 98;
		T[i] = 2;
		r[i] = 0.02;
		v[i] = 5;
	}

	for(int64_t i = 0; i < VECTOR_WIDTH; i++) {
		mask[i] = 1 % 2;
	}

	double sum = 0;	
	for(int n = 0; n < ROUNDS; n++) {
		for(int64_t i = 0; i < LONGVECTOR_LENGTH; i += VECTOR_WIDTH) {
			body(&S[i], &X[i], &T[i], &r[i], &v[i], &sum);
		}
	}

	std::cout << sum/(LONGVECTOR_LENGTH * ROUNDS) << std::endl;

	return 0;
}
