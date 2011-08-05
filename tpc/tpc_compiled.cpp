
#include <math.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <stdint.h>
#include <sys/time.h>
#include <iostream>

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

double* load(std::string name) {
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

void print6(std::string name, double* v) {
	std::cout << name << ": " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ", " << v[4] << ", " << v[5] << std::endl;
}

void print6(std::string name, int64_t* v) {
	std::cout << name << ": " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ", " << v[4] << ", " << v[5] << std::endl;
}


int main(int argc, char** argv) {

	int N = 6001215;

	double* l_quantity = load("quantity.bin");
	double* l_extendedprice = load("extendedprice.bin");
	double* l_discount = load("discount.bin");
	double* l_tax = load("tax.bin");
	
	int64_t* l_returnflag = load_int("returnflag.bin");
	int64_t* l_linestatus = load_int("linestatus.bin");

	double sq[6], sbp[6], sdp[6], sc[6], aq[6], ap[6], ad[6];
	int64_t count[6];

	for(int i = 0; i < 6; i++) {
		sq[i] = 0;
		sbp[i] = 0;
		sdp[i] = 0;
		sc[i] = 0;
		aq[i] = 0;
		ap[i] = 0;
		ad[i] = 0;
		count[i] = 0;
	}

	start_timing();

	for(int i = 0; i < N; i++) {
		int64_t group = (l_returnflag[i]+1) * (l_linestatus[i]+1) - 1;
		sq[group] += l_quantity[i];
		sbp[group] += l_extendedprice[i];
		sdp[group] += l_extendedprice[i] * (1-l_discount[i]);
		sc[group] += l_extendedprice[i] * (1-l_discount[i]) * (1+l_tax[i]);
		ad[group] += l_discount[i];
		count[group]++;
	}
	for(int i = 0; i < 6; i++) {
		aq[i] = sq[i]/count[i];
		ap[i] = sbp[i]/count[i];
		ad[i] = ad[i]/count[i];
	}

	print6("sum_quantity", sq);	
	print6("sum_base_price", sbp);	
	print6("sum_disc_price", sdp);	
	print6("sum_charge", sc);	
	print6("avg_quantity", aq);	
	print6("avg_price", ap);	
	print6("avg_disc", ad);
	print6("count", count);	

	printf("time: %f\n", end_timing());	
	
	return 0;
}
