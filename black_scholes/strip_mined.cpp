#include <stdio.h>
#include <math.h>

#define VW 256
#define VL 256
#define ROUNDS 20

template<int N>
void neg_op(double* i, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = -i[j];
	}
}

template<int N>
void abs_op(double* i, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = fabs(i[j]);
	}
}

template<int N>
void log_op(double* i, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = log(i[j]);
	}
}

template<int N>
void exp_op(double* i, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = exp(i[j]);
	}
}

template<int N>
void sqrt_op(double* i, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = sqrt(i[j]);
	}
}

template<int N>
void add_op(double* a, double* b, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = a[j] + b[j];
	}
}

template<int N>
void adds_op(double* a, double b, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = a[j] + b;
	}
}

template<int N>
void sub_op(double* a, double* b, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = a[j] - b[j];
	}
}

template<int N>
void mul_op(double* a, double* b, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = a[j] * b[j];
	}
}

template<int N>
void muls_op(double* a, double b, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = a[j] * b;
	}
}

template<int N>
void div_op(double* a, double* b, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = a[j] / b[j];
	}
}

template<int N>
void rep_op(double d, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = d;
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
	for(int j = 0; j < N; j++) {
		o[j] = i[j] > 0 ? 1.0 : 0.0;
	}
}

template<int N>
void sel_op(double* s, double* a, double* b, double* o) {
	for(int j = 0; j < N; j++) {
		o[j] = s[j] > 0 ? a[j] : b[j];
	}
}

template<int N>
double* getV() {
	return new double[N];
}

double* L = getV<VW>();
double* k = getV<VW>();
double* t0 = getV<VW>();
double* t1 = getV<VW>();
double* k2 = getV<VW>();
			double* k3 = getV<VW>();
			double* k4 = getV<VW>();
			double* k5 = getV<VW>();
			double* is2pi = getV<VW>();
			double* w = getV<VW>();

double* CND(double* X) {
			// CND
			abs_op<VW>(X, L);
			//rep_op<VW>(0.2316419, t0);
			//mul_op<VW>(L, t0, t0);
			muls_op<VW>(L, 0.2316419, t0);
			//rep_op<VW>(1.0, t1);
			//add_op<VW>(t0, t1, t0);
			adds_op<VW>(t0, 1.0, t0);
			rep_op<VW>(1.0, t1);
			div_op<VW>(t1, t0, k);

			mul_op<VW>(k, k, k2);
			mul_op<VW>(k2, k, k3);
			mul_op<VW>(k2, k2, k4);
			mul_op<VW>(k2, k3, k5);

			//rep_op<VW>(0.39894228040, is2pi);

			//rep_op<VW>(0.31938153, t0);
			//mul_op<VW>(t0, k, t0);
			muls_op<VW>(k, 0.31938153, t0);
			//rep_op<VW>(-0.356563782, t1);
			//mul_op<VW>(t1, k2, t1);
			muls_op<VW>(k, -0.356563782, t1);
			add_op<VW>(t0, t1, t0);
			//rep_op<VW>(1.781477937, t1);
			//mul_op<VW>(t1, k3, t1);
			muls_op<VW>(k, 1.781477937, t1);
			add_op<VW>(t0, t1, t0);
			//rep_op<VW>(-1.821255978, t1);
			//mul_op<VW>(t1, k4, t1);
			muls_op<VW>(k, -1.821255978, t1);
			add_op<VW>(t0, t1, t0);
			//rep_op<VW>(1.330274429, t1);
			//mul_op<VW>(t1, k5, t1);
			muls_op<VW>(k, 1.330274429, t1);
			add_op<VW>(t0, t1, w);
	
			//mul_op<VW>(w, is2pi, w);
			muls_op<VW>(w, 0.39894228040, w);
			neg_op<VW>(L, t0);
			mul_op<VW>(L, t0, t0);
			//rep_op<VW>(0.5, t1);
			//mul_op<VW>(t0, t1, t0);
			muls_op<VW>(t0, 0.5, t0);
			exp_op<VW>(t0, t0);
			mul_op<VW>(w, t0, w);

			gt0_op<VW>(X, t0);
			rep_op<VW>(1, t1);
			sub_op<VW>(t1, w, t1);
			sel_op<VW>(t0, t1, w, w);

	return w;
}

	double* S = getV<VW>();
	double* X = getV<VW>();
	double* T = getV<VW>();
	double* r = getV<VW>();
	double* v = getV<VW>();
	double* s0 = getV<VW>();
	double* s1 = getV<VW>();
	double* s2 = getV<VW>();
	double* d1 = getV<VW>();
	double* d2 = getV<VW>();
	double* cndd1 = CND(d1);
	double* result = getV<VW>();
	
double loop_body() {
	rep_op<VW>(100, S);
	rep_op<VW>(98, X);
	rep_op<VW>(2, T);
	rep_op<VW>(0.02, r);
	rep_op<VW>(5, v);

	div_op<VW>(S, X, s0);
	log_op<VW>(s0, s0);
	rep_op<VW>(log(10), s1);
	div_op<VW>(s0, s1, s0);
	

	mul_op<VW>(v, v, s1);
	//rep_op<VW>(0.5, s2);
	//mul_op<VW>(s1, s2, s1);
	muls_op<VW>(s1, 0.5, s1);

	add_op<VW>(s1, r, s1);
	mul_op<VW>(s1, T, s1);

	add_op<VW>(s0, s1, s0);
	
	sqrt_op<VW>(T, s1);
	mul_op<VW>(v, s1, s1);

	div_op<VW>(s0, s1, d1);

	sqrt_op<VW>(T, s0);
	mul_op<VW>(v, s0, s0);
	sub_op<VW>(d1, s0, d2);


	mul_op<VW>(S, CND(d1), s0);
	neg_op<VW>(r, s1);
	mul_op<VW>(s1, T, s1);
	exp_op<VW>(s1, s1);
	mul_op<VW>(X, s1, s1);
	mul_op<VW>(X, CND(d2), s1);
	sub_op<VW>(s0, s1, result);
	
	return sum_op<VW>(result);
}


int main(int argc, char** argv) {

	double acc = 0;
	for(int i = 0; i < ROUNDS; i++) {
		for(int j = 0; j < VL; j++)
			acc += loop_body();
	}

	printf("Result: %f\n", acc / (ROUNDS * VW * VL));

	return 0;
}

