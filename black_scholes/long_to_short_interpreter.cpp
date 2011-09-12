 #include <stdio.h>
#include <math.h>
#include <vector>
#include<string>
#include<sstream>

#define VL (LENGTH)


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
_(sum,__) \
_(gt0,u) \
_(sel,__) \
/* for trace interpreter only */ \
_(ld,u)


#define MAKE_ENUM(op,t) bc_##op,
enum ByteCode {
	INTERPRETER_OPS(MAKE_ENUM)
	INST_SIZE
};
#undef MAKE_ENUM

const char * bytecode_to_string(ByteCode bc) {
	#define MAKE_STRINGS(op,t) #op ,
	static const char * strs[] = { INTERPRETER_OPS(MAKE_STRINGS) 0 };
	#undef MAKE_STRINGS
	return strs[bc];
}

enum InstType {
	I_u,
	I_b,
	I_bs,
	I___
};
InstType bytecode_to_type(ByteCode bc) {
	#define MAKE_TYPE(op,t) I_##t,
	static InstType strs[] = { INTERPRETER_OPS(MAKE_TYPE) I___ };
	#undef MAKE_TYPE
	return strs[bc];
}

struct Inst {
	ByteCode code;
	int r, a, b;
};

struct Vector {
	int ref_count;
	int size;
	double data[0];
	std::string toString() {
		std::ostringstream out;
		out << "[" << size << "] {";
		for(int i = 0; i < 10; i++)
			out << data[i] << ((i + 1 == size) ? "}" : ",");
		if(size > 10)
			out << " ...";
		return out.str();
	}
};
enum ValueType {
	T_NULL, /* uninitialized value */
	T_VECTOR,
	T_FUTURE
};

struct Value {
	union {
		Vector * v;
		uint64_t idx;
	};
	ValueType typ;
	std::string toString() {
		std::ostringstream out;
		if(typ == T_NULL)
			out << "null";
		else if(typ == T_FUTURE)
			out << "f" << idx;
		else {
			out << v->toString();
		}
		return out.str();
	}
};

Value registers[32];

std::vector<double> constants;
std::vector<Inst> code;

struct IRNode {
	ByteCode bc;
	int a,b;
	int output; // >0 if this IRNode an output of the trace, indicating where to write the result
	double * reg;
};

#define MAX_NODES 128
#define MAX_INPUTS 32

struct Trace {
	Vector * inputs[MAX_INPUTS];
	int input_nodes[MAX_INPUTS];
	
	IRNode nodes[MAX_NODES];
	int n_nodes;
	int n_inputs;
	int max_size;
	int iteration;
	double registers[MAX_NODES][BLOCK];
};
Trace trace;


void Vector_acquire(Vector * v) {
	v->ref_count++;
}
void Vector_release(Vector * v) {
	v->ref_count--;
	if(v->ref_count == 0)
		free(v);
}

void Vector_create(Vector ** v, int size) {
	*v = (Vector *) malloc(sizeof(Vector) + sizeof(double) * size);
	**v = (Vector) { 1, size};
}


void Value_assign(Value * lhs, Vector * rhs) {
	if(lhs->typ == T_VECTOR) {
		Vector_release(lhs->v);
	} else if(lhs->typ == T_FUTURE) {
		trace.nodes[lhs->idx].output = -1;
	}
	lhs->typ = T_VECTOR;
	lhs->v = rhs;
}
void Value_assign(Value * lhs, int future) {
	if(lhs->typ == T_VECTOR) {
		Vector_release(lhs->v);
	} else if(lhs->typ == T_FUTURE) {
		trace.nodes[lhs->idx].output = -1;
	}
	lhs->typ = T_FUTURE;
	lhs->idx = future;
}



int insert_const(double d) {
	int r = constants.size();
	constants.push_back(d);
	return r;
}
void insert_inst(ByteCode c, int r, int a, int b) {
	Inst i = {c,r,a,b};
	code.push_back(i);
}

void trace_flush();

void trace_reserve(int inst, int inputs) {
	if(trace.n_nodes + inst > MAX_NODES || trace.n_inputs + inputs > MAX_INPUTS)
		trace_flush();
}
int trace_record(ByteCode bc, int r, int a, int b) {
	IRNode & n = trace.nodes[trace.n_nodes];
	n.bc = bc;
	n.a = a;
	n.b = b;
	n.output = r;
	n.reg = trace.registers[trace.n_nodes]; //for now each node gets its own temporary space, we should really register allocate this to keep the working set as small as possible
	return trace.n_nodes++;
}

int trace_input(Vector * v) {
	//inputs are often repeated locally
	//to avoid redundant loads we look back in the input list a bit to see if we can find one
	//an alternative would be to modify all Values loaded into the interpreter to have type T_FUTURE and point to the ld node already in the interpreter
	int look_back = std::min(trace.n_inputs,3);
	for(int i = trace.n_inputs - look_back; i < trace.n_inputs; i++)
		if(v == trace.inputs[i])
			return trace.input_nodes[i];
	Vector_acquire(v);
	trace.max_size = std::max(trace.max_size,v->size);
	trace.inputs[trace.n_inputs] = v;
	return trace.input_nodes[trace.n_inputs] = trace_record(bc_ld,-1,trace.n_inputs++,0);
}


void trace_interp();

void trace_flush() {
#ifdef DEBUG_INTERP
	printf("inputs [%d]\n",trace.n_inputs);
	for(int i = 0; i < trace.n_inputs; i++) {
		printf("i%d = %s\n",i,trace.inputs[i]->toString().c_str());
	}
	printf("insts [%d]\n",trace.n_nodes);
	for(int i = 0; i < trace.n_nodes; i++) {
		IRNode & node = trace.nodes[i];
		InstType it = bytecode_to_type(node.bc);
		if(it == I_u)
			printf("[%d] r%d = %s r%d\n",node.output,i,bytecode_to_string(node.bc),node.a);
		else if(it == I_b)
			printf("[%d] r%d = %s r%d r%d\n",node.output,i,bytecode_to_string(node.bc),node.a,node.b);
		else if(it == I_bs)
			printf("[%d] r%d = %s r%d %f\n",node.output,i,bytecode_to_string(node.bc),node.a,constants[node.b]);
	}
#endif
	trace_interp();
}
#define MAKE_u(op) void op##_ins(int i, int o) { /*printf("%s %d %d\n", #op, i, o);*/ insert_inst(bc_##op,o,i,0); }
#define MAKE_b(op) void op##_ins(int a, int b, int o) { /*printf("%s %d %d %d\n", #op, a, b, o);*/ insert_inst(bc_##op,o,a,b);}
#define MAKE_bs(op) void op##_ins(int a, double b, int o) { /*printf("%s %d %f %d\n", #op, a, b, o);*/ insert_inst(bc_##op,o,a,insert_const(b));}
#define MAKE___(op)

#define MAKE_INS(op,typ) MAKE_##typ(op)
INTERPRETER_OPS(MAKE_INS)

void reps_ins(double a, double b, int o) { insert_inst(bc_reps,o,insert_const(a),insert_const(b));}
void sum_ins(int i, int o) { insert_inst(bc_sum,o,i,0); }
void sel_ins(int a, int b, int o) { insert_inst(bc_sel,o,a,b);}

#undef MAKE_u
#undef MAKE_b
#undef MAKE_bs
#undef MAKE___
#undef MAKE_INS

#ifdef DEBUG_PRINT
void PR2(const char * op, double * value) {
	printf("%s %f\n",op,value[0]);
}
#else
#define PR2(a,b) do{}while(0);
#endif

#define UNARY_OP(op, impl) \
void op##_op(Vector* i, Value * optr) { \
	Vector * o; \
	Vector_create(&o,i->size); \
	for(int j = 0; j < i->size; j++) { \
		double a = i->data[j];\
		o->data[j] = impl;\
	} \
	Value_assign(optr,o); \
} \
void op##_top(IRNode * node) { \
	double * o = node->reg; \
	double * i = trace.nodes[node->a].reg; \
	for(int j = 0; j < BLOCK; j++) { \
		double a = i[j]; \
		o[j] = impl; \
	} \
	PR2(#op,o);\
}

UNARY_OP(neg,-a)
UNARY_OP(abs,fabs(a))
UNARY_OP(recip,1.0 / a)
UNARY_OP(log,log(a))
UNARY_OP(exp,exp(a))
UNARY_OP(sqrt,sqrt(a))
UNARY_OP(gt0,a > 0 ? 1.0 : 0.0)

void ld_op(Vector* i, Value * optr)  {}
void ld_top(IRNode * node) {
	Vector * v = trace.inputs[node->a];
	node->reg = v->data + trace.iteration;
}

#define BINARY_OP(op,impl) \
void op##_op(Vector* av, Vector* bv, Value * optr) { \
	Vector * o; \
	Vector_create(&o,av->size); \
	/*NYI: vectors of different size*/\
	for(int j = 0; j < av->size; j++) { \
		double a = av->data[j];\
		double b = bv->data[j];\
		o->data[j] = impl;\
		impl;\
	} \
	Value_assign(optr,o); \
}\
void op##_top(IRNode * node) { \
	double * o = node->reg; \
	double * av = trace.nodes[node->a].reg; \
	double * bv = trace.nodes[node->b].reg; \
	for(int j = 0; j < BLOCK; j++) { \
		double a = av[j]; \
		double b = bv[j]; \
		o[j] = impl; \
	} \
	PR2(#op,o);\
}

BINARY_OP(add,a + b)
BINARY_OP(sub,a - b)
BINARY_OP(mul,a * b)
BINARY_OP(div,a / b)


#define BINARY_SCALAR_OP(op,impl) \
void op##_op(Vector* av, double b, Value * optr) { \
	Vector * o; \
	Vector_create(&o,av->size); \
	/*NYI: vectors of different size*/\
	for(int j = 0; j < av->size; j++) { \
		double a = av->data[j];\
		o->data[j] = impl;\
	} \
	Value_assign(optr,o); \
}\
void op##_top(IRNode * node) { \
	double * o = node->reg; \
	double * av = trace.nodes[node->a].reg; \
	double b = constants[node->b]; \
	for(int j = 0; j < BLOCK; j++) { \
		double a = av[j]; \
		o[j] = impl; \
	} \
	PR2(#op,o);\
}

BINARY_SCALAR_OP(adds,a + b)
BINARY_SCALAR_OP(muls,a * b)
BINARY_SCALAR_OP(divs,a / b)


void reps_op(double v, double size, Value * optr) {
	Vector * o;
	Vector_create(&o,size);
	for(int j = 0; j < o->size; j++) {
		o->data[j] = v;
	}
	Value_assign(optr,o);
}
void reps_top(IRNode *) {}
void sum_top(IRNode *) {}
void sel_top(IRNode *) {}

void sum_op(Vector* i, Value * optr) {
	Vector * o;
	Vector_create(&o,1); 
	o->data[0] = 0.0;
	for(int j = 0; j < i->size; j++) {
		o->data[0] += i->data[j];
	}
	Value_assign(optr,o);
}

void sel_op(Vector* s, Vector* a, Vector* b) {
	for(int j = 0; j < b->size; j++) {
		b->data[j] = s->data[j] > 0 ? a->data[j] : b->data[j];
	}
}


//#define DEBUG_INTERP

#ifdef DEBUG_INTERP
void PR(const char * op, int i) {
	Value & rhs = registers[i];
	printf("%s: %s\n",op,rhs.toString().c_str());
}
#else
#define PR(a,b) do{} while(0)
#endif


#define MAKE_u(op) \
int op##_interp(Inst * i) {\
	Value & r = registers[i->r];\
	Value & a = registers[i->a];\
	if(a.typ == T_FUTURE) { \
		trace_reserve(1,0); \
		int v = trace_record(bc_##op,i->r,a.idx,0); \
		Value_assign(&r,v); \
	} else if(a.v->size > BLOCK) { \
		trace_reserve(2,1); \
		int idx = trace_input(a.v); \
		int v = trace_record(bc_##op,i->r,idx,0); \
		Value_assign(&r,v);  \
	} else \
		op##_op(a.v,&r);\
	PR(#op,i->r);\
	return 1;\
}

#define MAKE_b(op) int op##_interp(Inst * i) { \
	Value & r = registers[i->r];\
	Value & a = registers[i->a];\
	Value & b = registers[i->b];\
	if(a.typ == T_FUTURE) {\
		if(b.typ == T_FUTURE) {\
			trace_reserve(1,0);\
			Value_assign(&r,trace_record(bc_##op,i->r,a.idx,b.idx));\
		} else {\
			trace_reserve(2,1);\
			int idx = trace_input(b.v);\
			Value_assign(&r,trace_record(bc_##op,i->r,a.idx,idx));\
		}\
	} else if(a.v->size > BLOCK) {\
		if(b.typ == T_FUTURE) {\
			trace_reserve(2,1);\
			int idx = trace_input(a.v);\
			Value_assign(&r,trace_record(bc_##op,i->r,idx,b.idx));\
		} else {\
			trace_reserve(3,2);\
			int aidx = trace_input(a.v);\
			int bidx = trace_input(b.v);\
			Value_assign(&r,trace_record(bc_##op,i->r,aidx,bidx));\
		}\
	} else {\
		if(b.typ == T_FUTURE) {\
			trace_reserve(2,1);\
			int idx = trace_input(a.v);\
			Value_assign(&r,trace_record(bc_##op,i->r,idx,b.idx));\
		} else if (b.v->size > BLOCK) {\
			trace_reserve(3,2);\
			int aidx = trace_input(a.v);\
			int bidx = trace_input(b.v);\
			Value_assign(&r,trace_record(bc_##op,i->r,aidx,bidx));\
		} else {\
			op##_op(a.v,b.v,&r);\
		} \
	} \
	PR(#op,i->r); \
	return 1;\
}
#define MAKE_bs(op) int op##_interp(Inst * i) { \
	Value & r = registers[i->r];\
	Value & a = registers[i->a];\
	if(a.typ == T_FUTURE) { \
		trace_reserve(1,0); \
		int v = trace_record(bc_##op,i->r,a.idx,i->b); \
		Value_assign(&r,v); \
	} else if(a.v->size > BLOCK) { \
		trace_reserve(2,1); \
		int idx = trace_input(a.v); \
		int v = trace_record(bc_##op,i->r,idx,i->b); \
		Value_assign(&r,v);  \
	} else \
		op##_op(a.v,constants[i->b],&r);\
	PR(#op,i->r); \
	return 1;\
}
#define MAKE___(op)

#define MAKE_INTERP(op,typ) MAKE_##typ(op)
INTERPRETER_OPS(MAKE_INTERP)
int reps_interp(Inst * i) { reps_op(constants[i->a],constants[i->b],&registers[i->r]); PR("reps",i->r); return 1;}
int sum_interp(Inst * i) {
	Value & v = registers[i->a];
	if(v.typ == T_FUTURE)
		trace_flush();
	sum_op(v.v,&registers[i->r]);
	PR("sum",i->r);
	return 1;
}
int sel_interp(Inst * i) {
	Value & r = registers[i->r];
	Value & a = registers[i->a];
	Value & b = registers[i->b];
	if(r.typ == T_FUTURE || a.typ == T_FUTURE || b.typ == T_FUTURE)
		trace_flush();
	sel_op(a.v,b.v,r.v);
	PR("sel",i->r);
	return 1;
}

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
		default: ip -= 1; break; //make sure it can't unroll this loop
		}
	}
	
}
void trace_interp() {
	static int trace_count = 0;
	int outputs_list[][4] = { { 9, 2, 3,19}, {9, 2 ,3,15}, {21,0,0,0} }; //we need a way to determine outputs better but for now I am just hand-coding them
	int outputs_lengths[] = { 4, 4, 1 };

	int * outputs = outputs_list[trace_count];
	int n_outputs = outputs_lengths[trace_count];
	trace_count = (trace_count + 1) % 3;

	Vector * output_values[4];


	//allocate output vectors, allocate vectors
	for(int i = 0; i < n_outputs; i++) {
		Vector_create(&output_values[i],trace.max_size);
		int idx = registers[outputs[i]].idx;
		trace.nodes[idx].reg = output_values[i]->data;
	}

	for(trace.iteration = 0; trace.iteration < trace.max_size; trace.iteration += BLOCK) {

		for(int i = 0; i < trace.n_nodes; i++) {
			IRNode & node = trace.nodes[i];
			switch(node.bc) {
#define INTERP_CASE(op,_) case bc_##op: op##_top(&node); break;
		INTERPRETER_OPS(INTERP_CASE)
#undef INTERP_CASE
			}
		}
		//advance the output temporaries
		for(int i = 0; i < n_outputs; i++) {
			int idx = registers[outputs[i]].idx;
			trace.nodes[idx].reg += BLOCK;
		}
	}

	for(int i = 0; i < n_outputs; i++) {
		Value_assign(&registers[outputs[i]],output_values[i]);
	}
	for(int i = 0; i < trace.n_inputs; i++)
		Vector_release(trace.inputs[i]);

	//reset the trace variables
	trace.max_size = 0;
	trace.n_inputs = 0;
	trace.n_nodes = 0;
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
		acc += registers[ret_val].v->data[0];
	}

	printf("Result: %f\n", acc / (ROUNDS * VL));

	return 0;
}

