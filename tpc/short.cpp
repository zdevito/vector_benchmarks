
#include <iostream>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sys/time.h>

#define VECTOR_WIDTH 4
#define NOINLINE __attribute__ ((noinline))

uint64_t readTime()
{
  timeval time_tt;
  gettimeofday(&time_tt, NULL);
  return (uint64_t)time_tt.tv_sec * 1000 * 1000 + (uint64_t)time_tt.tv_usec;
}


static uint64_t start_time;

void start_timing() {
	start_time = readTime();
}
double end_timing() {
	uint64_t end_time = readTime();
	return (end_time - start_time) / (1000.0 * 1000.0);
}

double* load_double(std::string name) {
	std::ifstream file(name.c_str(), std::ios::in | std::ios::binary);
	double* result = new double[6001215];
	file.read((char*)&result[0], sizeof(double) * 6001215);
	file.close();
	return result;
}

int64_t* load_int(std::string name) {
	std::ifstream file(name.c_str(), std::ios::in | std::ios::binary);
	int64_t* result = new int64_t[6001215];
	file.read((char*)&result[0], sizeof(int64_t) * 6001215);
	file.close();
	return result;
}

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
NOINLINE void adds_d(double*__restrict__ a, double*__restrict__ d, double b) { APPLY(a[i] = d[i]+b) }

NOINLINE void sub_d(double*__restrict__ a, double*__restrict__ b) { APPLY(a[i] -= b[i]) }
NOINLINE void subs_d(double*__restrict__ a, double d, double*__restrict__ b) { APPLY(a[i] = d-b[i]) }
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

NOINLINE void copy_i(int64_t*__restrict__ a, int64_t*__restrict__ b) { APPLY(a[i] = b[i]) }
NOINLINE void adds_i(int64_t*__restrict__ a, int64_t b) { APPLY(a[i] += b) }
NOINLINE void adds_i(int64_t*__restrict__ a, int64_t*__restrict__ d, int64_t b) { APPLY(a[i] = d[i] + b) }
NOINLINE void subs_i(int64_t*__restrict__ a, int64_t b) { APPLY(a[i] -= b) }
NOINLINE void mul_i(int64_t*__restrict__ a, int64_t*__restrict__ b) { APPLY(a[i] = a[i]*b[i]) }
NOINLINE void muls_i(int64_t*__restrict__ a, int64_t b) { APPLY(a[i] *= b) }

NOINLINE void count_i(int64_t*__restrict__ a, int64_t*__restrict__ group) { APPLY(a[group[i]]++) }

NOINLINE void cast_id(double*__restrict__ a, int64_t*__restrict__ b) { APPLY(a[i] = b[i]) }

double  d_reg[32][VECTOR_WIDTH] __attribute__((aligned (16)));
int64_t i_reg[32][VECTOR_WIDTH];

void print6(std::string name, double* v) {
	std::cout << name << ": " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ", " << v[4] << ", " << v[5] << std::endl;
}

void print6(std::string name, int64_t* v) {
	std::cout << name << ": " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ", " << v[4] << ", " << v[5] << std::endl;
}


int main(int argc, char** argv) {

	int64_t N = 6001215/VECTOR_WIDTH * VECTOR_WIDTH;  // have to handle non-multiple lengths

	double* l_quantity = load_double("quantity.bin");
	double* l_extendedprice = load_double("extendedprice.bin");
	double* l_discount = load_double("discount.bin");
	double* l_tax = load_double("tax.bin");
	
	int64_t* l_returnflag = load_int("returnflag.bin");
	int64_t* l_linestatus = load_int("linestatus.bin");

	int64_t* group = new int64_t[N];
	for(int i = 0; i < N; i++) {group[i] = (l_returnflag[i]+1) * (l_linestatus[i]+1) - 1;}

	// the output vectors, this would have to be allocated by the long-to-short
	double sq[VECTOR_WIDTH], sbp[VECTOR_WIDTH], sdp[VECTOR_WIDTH], sc[VECTOR_WIDTH], aq[VECTOR_WIDTH], ap[VECTOR_WIDTH], ad[VECTOR_WIDTH];
	int64_t cnt[VECTOR_WIDTH];

	for(int i = 0; i < VECTOR_WIDTH; i++) {
		sq[i] = 0;
		sbp[i] = 0;
		sdp[i] = 0;
		sc[i] = 0;
		aq[i] = 0;
		ap[i] = 0;
		ad[i] = 0;
		cnt[i] = 0;
	}

	int64_t *i0 = i_reg[0], *i1 = i_reg[1];
	double *d0 = d_reg[0], *d1 = d_reg[1];

	start_timing();

	for(int64_t i = 0; i < N; i+= VECTOR_WIDTH) {
		// compute groups first
		//copy_i(i0, &l_returnflag[i]);
		//adds_i(i0, 1);
		//copy_i(i1, &l_linestatus[i]);
		//adds_i(i1, 1);
		//mul_i(i0, i1);
		//subs_i(i0, 1);
		//copy_i(i0, &group[i]);
		adds_i(i0, &l_returnflag[i], 1);
		adds_i(i1, &l_linestatus[i], 1);
		mul_i(i0, i1);
		subs_i(i0, 1);

		// now, compute aggregates
		sum_d(sq, &l_quantity[i], i0);
		sum_d(sbp, &l_extendedprice[i], i0);
		//copy_d(d0, &l_discount[i]);
		//rsubs_d(d0, 1);
		subs_d(d0, 1, &l_discount[i]);
		mul_d(d0, &l_extendedprice[i]);
		sum_d(sdp, d0, i0);
		
		//copy_d(d1, &l_tax[i]);
		//adds_d(d1, 1);
		adds_d(d1, &l_tax[i], 1);
		mul_d(d0, d1);
		sum_d(sc, d0, i0);

		// computing average needs to be delayed until very end?
		sum_d(ad, &l_discount[i], i0);

		count_i(cnt, i0);
	}
	copy_d(aq, sq);
	copy_d(ap, sbp);
	cast_id(d0, cnt);
	div_d(aq, d0);
	div_d(ap, d0);
	div_d(ad, d0);

	print6("sum_quantity", sq);	
	print6("sum_base_price", sbp);	
	print6("sum_disc_price", sdp);	
	print6("sum_charge", sc);	
	print6("avg_quantity", aq);	
	print6("avg_price", ap);	
	print6("avg_disc", ad);
	print6("count", cnt);	

	printf("time: %f\n", end_timing());	

	return 0;
}

