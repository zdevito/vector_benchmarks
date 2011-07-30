 #include <stdio.h>
#include <math.h>
#include <vector>

#define VL (256 * 256)
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
_(divs,bs) \
_(recip,u) \
_(reps,__) \
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

struct Value {
	int ref_count;
	int size;
	double data[0];
};

void Value_acquire(Value * v) {
	v->ref_count++;
}
void Value_release(Value * v) {
	v->ref_count--;
	if(v->ref_count == 0)
		free(v);
}

void Value_create(Value ** v, int size) {
	*v = (Value *) malloc(sizeof(Value) + sizeof(double) * size);
	**v = (Value) { 1, size};
}

Value * registers[32];
void Value_assign(Value ** lhs, Value * rhs) {
	if(*lhs != NULL)
		Value_release(*lhs);
	*lhs = rhs;
}


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
#define MAKE___(op)

#define MAKE_INS(op,typ) MAKE_##typ(op)
INTERPRETER_OPS(MAKE_INS)

void reps_ins(double a, double b, int o) { insert_inst(bc_reps,o,insert_const(a),insert_const(b));}

#undef MAKE_u
#undef MAKE_b
#undef MAKE_us
#undef MAKE_bs
#undef MAKE___
#undef MAKE_INS




#define UNARY_OP(op, impl) \
void op##_op(Value* i, Value** optr) { \
	Value * o; \
	Value_create(&o,i->size); \
	for(int j = 0; j < i->size; j++) { \
		impl;\
	} \
	Value_assign(optr,o); \
}

UNARY_OP(neg,o->data[j] = -i->data[j])
UNARY_OP(abs,o->data[j] = fabs(i->data[j]))
UNARY_OP(recip,o->data[j] = 1.0 / i->data[j])
UNARY_OP(log,o->data[j] = log(i->data[j]))
UNARY_OP(exp,o->data[j] = exp(i->data[j]))
UNARY_OP(sqrt,o->data[j] = sqrt(i->data[j]))
UNARY_OP(gt0,o->data[j] = i->data[j] > 0 ? 1.0 : 0.0)

#define BINARY_OP(op,impl) \
void op##_op(Value* a, Value* b, Value** optr) { \
	Value * o; \
	Value_create(&o,a->size); \
	/*NYI: vectors of different size*/\
	for(int j = 0; j < a->size; j++) { \
		impl;\
	} \
	Value_assign(optr,o); \
}

BINARY_OP(add,o->data[j] = a->data[j] + b->data[j])
BINARY_OP(sub,o->data[j] = a->data[j] - b->data[j])
BINARY_OP(mul,o->data[j] = a->data[j] * b->data[j])
BINARY_OP(div,o->data[j] = a->data[j] / b->data[j])


#define BINARY_SCALAR_OP(op,impl) \
void op##_op(Value* a, double b, Value** optr) { \
	Value * o; \
	Value_create(&o,a->size); \
	/*NYI: vectors of different size*/\
	for(int j = 0; j < a->size; j++) { \
		impl;\
	} \
	Value_assign(optr,o); \
}

BINARY_SCALAR_OP(adds,o->data[j] = a->data[j] + b)
BINARY_SCALAR_OP(muls,o->data[j] = a->data[j] * b)
BINARY_SCALAR_OP(divs,o->data[j] = a->data[j] / b)


void reps_op(double v, double size, Value** optr) {
	Value * o;
	Value_create(&o,size);
	for(int j = 0; j < o->size; j++) {
		o->data[j] = v;
	}
	Value_assign(optr,o);
}

void sum_op(Value* i, Value** optr) {
	Value * o;
	Value_create(&o,1); 
	o->data[0] = 0.0;
	for(int j = 0; j < i->size; j++) {
		o->data[0] += i->data[j];
	}
	Value_assign(optr,o);
}

void sel_op(Value* s, Value* a, Value** bptr) {
	Value * b = *bptr;
	for(int j = 0; j < b->size; j++) {
		b->data[j] = s->data[j] > 0 ? a->data[j] : b->data[j];
	}
}



//#define DEBUG_INTERP

#ifdef DEBUG_INTERP
void PR(const char * op, int i) {
	Value * rhs = registers[i];
	printf("%s v[%d][%d] = { ",op,i,rhs->size);
	for(int j = 0; j < rhs->size; j++)
		printf("%f ",rhs->data[j]);
	printf("}\n");
}
#else
#define PR(a,b) do{} while(0)
#endif

#define MAKE_u(op) int op##_interp(Inst * i) { op##_op(registers[i->a],&registers[i->r]); PR(#op,i->r); return 1;}
#define MAKE_b(op) int op##_interp(Inst * i) { op##_op(registers[i->a],registers[i->b],&registers[i->r]); PR(#op,i->r); return 1;}
#define MAKE_us(op) int op##_interp(Inst * i) { op##_op(constants[i->a],&registers[i->r]); PR(#op,i->r); return 1;}
#define MAKE_bs(op) int op##_interp(Inst * i) { op##_op(registers[i->a],constants[i->b],&registers[i->r]); PR(#op,i->r); return 1; }
#define MAKE___(op)

#define MAKE_INTERP(op,typ) MAKE_##typ(op)
INTERPRETER_OPS(MAKE_INTERP)
int reps_interp(Inst * i) { reps_op(constants[i->a],constants[i->b],&registers[i->r]); PR("reps",i->r); return 1;}

#undef MAKE_u
#undef MAKE_b
#undef MAKE_us
#undef MAKE_bs
#undef MAKE___
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
			recip_ins(t0,k);

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
			neg_ins(w,t1);
			adds_ins(t1, 1, t1);
			sel_ins(t0, t1, w);

	return w;
}
	
void loop_body() {
	reps_ins(100, VL, S);
	reps_ins(98, VL, X);
	reps_ins(2, VL, T);
	reps_ins(0.02,VL, r);
	reps_ins(5,VL, v);

	div_ins(S, X, s0);
	log_ins(s0, s0);
	divs_ins(s0,log(10), s0);

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
	loop_body();
	
	double acc = 0;
	for(int i = 0; i < ROUNDS; i++) {
		interp();
		acc += registers[ret_val]->data[0];
	}

	printf("Result: %f\n", acc / (ROUNDS * VL));

	return 0;
}

