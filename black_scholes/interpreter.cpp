 #include <stdio.h>
#include <math.h>
#include <vector>

#define VW 256
#define VL 256
#define ROUNDS 20


#define INTERPRETER_OPS(_) \
_(neg,u) \
_(abs,u) \
_(log,u) \
_(exp,u) \
_(sqrt,u) \
_(add,b) \
_(adds,bs) \
_(sub,b) \
_(mul,b) \
_(muls,bs) \
_(div,b) \
_(rep,us) \
_(sum,u) \
_(gt0,u) \
_(sel,b) \


#define MAKE_ENUM(op,t) bc_##op,
enum ByteCode {
	INTERPRETER_OPS(MAKE_ENUM)
	INST_SIZE
};
#undef MAKE_ENUM

struct Inst {
	ByteCode code;
	int r, a, b;
};

double * registers[32];

std::vector<double> constants;
std::vector<Inst> code;

int insert_const(double d) {
	int r = constants.size();
	constants.push_back(d);
	return r;
}
void insert_inst(ByteCode c, int r, int a, int b) {
	Inst i = {c,r,a,b};
	code.push_back(i);
}

#define MAKE_u(op) void op##_ins(int i, int o) { /*printf("%s %d %d\n", #op, i, o);*/ insert_inst(bc_##op,o,i,0); }
#define MAKE_b(op) void op##_ins(int a, int b, int o) { /*printf("%s %d %d %d\n", #op, a, b, o);*/ insert_inst(bc_##op,o,a,b);}
#define MAKE_us(op) void op##_ins(double i, int o) { /*printf("%s %f %d\n", #op, i, o);*/ insert_inst(bc_##op,o,insert_const(i),0);}
#define MAKE_bs(op) void op##_ins(int a, double b, int o) { /*printf("%s %d %f %d\n", #op, a, b, o);*/ insert_inst(bc_##op,o,a,insert_const(b));}

#define MAKE_INS(op,typ) MAKE_##typ(op)
INTERPRETER_OPS(MAKE_INS)

#undef MAKE_u
#undef MAKE_b
#undef MAKE_us
#undef MAKE_bs
#undef MAKE_INS



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
void sum_op(double* i, double * o) {
	double s = 0;
	for(int j = 0; j < N; j++) {
		s += i[j];
	}
	o[0] = s;
}

template<int N>
void gt0_op(double* i, double*o ) {
	for(int j = 0; j < N; j++) {
		o[j] = i[j] > 0 ? 1.0 : 0.0;
	}
}

template<int N>
void sel_op(double* s, double* a, double* b) {
	for(int j = 0; j < N; j++) {
		b[j] = s[j] > 0 ? a[j] : b[j];
	}
}


#define MAKE_u(op) int op##_interp(Inst * i) { op##_op<VW>(registers[i->a],registers[i->r]); return 1;}
#define MAKE_b(op) int op##_interp(Inst * i) { op##_op<VW>(registers[i->a],registers[i->b],registers[i->r]); return 1;}
#define MAKE_us(op) int op##_interp(Inst * i) { op##_op<VW>(constants[i->a],registers[i->r]); return 1;}
#define MAKE_bs(op) int op##_interp(Inst * i) { op##_op<VW>(registers[i->a],constants[i->b],registers[i->r]); return 1; }

#define MAKE_INTERP(op,typ) MAKE_##typ(op)
INTERPRETER_OPS(MAKE_INTERP)

#undef MAKE_u
#undef MAKE_b
#undef MAKE_us
#undef MAKE_bs
#undef MAKE_INS

void interp() {
	Inst * ip = &code[0];
	Inst * end = ip + code.size();
	while(ip != end) {
		switch(ip->code) {
#define INTERP_CASE(op,_) case bc_##op: ip += op##_interp(ip); break;
		INTERPRETER_OPS(INTERP_CASE)
#undef INTERP_CASE
		default: ip -= 1; //make sure it can't unroll this loop
		}
	}
	
}

template<int N>
double* getV() {
	return new double[N];
}
enum {
	L,k,t0,t1,k2,k3,k4,k5,is2pi,w,
	S,X,T,r,v,s0,s1,s2,d1,d2,cndd1,result,ret_val
};

int CND(int X) {
			// CND
			abs_ins(X, L);
			//rep_ins(0.2316419, t0);
			//mul_ins(L, t0, t0);
			muls_ins(L, 0.2316419, t0);
			//rep_ins(1.0, t1);
			//add_ins(t0, t1, t0);
			adds_ins(t0, 1.0, t0);
			rep_ins(1.0, t1);
			div_ins(t1, t0, k);

			mul_ins(k, k, k2);
			mul_ins(k2, k, k3);
			mul_ins(k2, k2, k4);
			mul_ins(k2, k3, k5);

			//rep_ins(0.39894228040, is2pi);

			//rep_ins(0.31938153, t0);
			//mul_ins(t0, k, t0);
			muls_ins(k, 0.31938153, t0);
			//rep_ins(-0.356563782, t1);
			//mul_ins(t1, k2, t1);
			muls_ins(k, -0.356563782, t1);
			add_ins(t0, t1, t0);
			//rep_ins(1.781477937, t1);
			//mul_ins(t1, k3, t1);
			muls_ins(k, 1.781477937, t1);
			add_ins(t0, t1, t0);
			//rep_ins(-1.821255978, t1);
			//mul_ins(t1, k4, t1);
			muls_ins(k, -1.821255978, t1);
			add_ins(t0, t1, t0);
			//rep_ins(1.330274429, t1);
			//mul_ins(t1, k5, t1);
			muls_ins(k, 1.330274429, t1);
			add_ins(t0, t1, w);
	
			//mul_ins(w, is2pi, w);
			muls_ins(w, 0.39894228040, w);
			neg_ins(L, t0);
			mul_ins(L, t0, t0);
			//rep_ins(0.5, t1);
			//mul_ins(t0, t1, t0);
			muls_ins(t0, 0.5, t0);
			exp_ins(t0, t0);
			mul_ins(w, t0, w);

			gt0_ins(X, t0);
			rep_ins(1, t1);
			sub_ins(t1, w, t1);
			sel_ins(t0, t1, w);

	return w;
}
	
void loop_body() {
	rep_ins(100, S);
	rep_ins(98, X);
	rep_ins(2, T);
	rep_ins(0.02, r);
	rep_ins(5, v);

	div_ins(S, X, s0);
	log_ins(s0, s0);
	rep_ins(log(10), s1);
	div_ins(s0, s1, s0);
	

	mul_ins(v, v, s1);
	//rep_ins(0.5, s2);
	//mul_ins(s1, s2, s1);
	muls_ins(s1, 0.5, s1);

	add_ins(s1, r, s1);
	mul_ins(s1, T, s1);

	add_ins(s0, s1, s0);
	
	sqrt_ins(T, s1);
	mul_ins(v, s1, s1);

	div_ins(s0, s1, d1);

	sqrt_ins(T, s0);
	mul_ins(v, s0, s0);
	sub_ins(d1, s0, d2);


	mul_ins(S, CND(d1), s0);
	neg_ins(r, s1);
	mul_ins(s1, T, s1);
	exp_ins(s1, s1);
	mul_ins(X, s1, s1);
	mul_ins(X, CND(d2), s1);
	sub_ins(s0, s1, result);
	
	sum_ins(result,ret_val);
}

int main(int argc, char** argv) {

	for(int i = 0; i < 32; i++) {
		registers[i] = getV<VW>();
	}
	
	loop_body();
	
	double acc = 0;
	for(int i = 0; i < ROUNDS; i++) {
		for(int j = 0; j < VL; j++) {
			interp();
			acc += registers[ret_val][0];
		}
	}

	printf("Result: %f\n", acc / (ROUNDS * VW * VL));

	return 0;
}

