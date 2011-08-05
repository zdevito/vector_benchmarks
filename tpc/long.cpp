
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

double sum_op(int N, double* i) {
	double s = 0;
	for(int j = 0; j < N; j++) {
		s += i[j];
	}
	return s;
}

double* sum_op(int N, double* i, int64_t* group) {
	double* result = new double[6];
	for(int j = 0; j < 6; j++) result[j] = 0;
	for(int j = 0; j < N; j++) {
		result[group[j]] += i[j];
	}
	return result;
}

double* count_op(int N, int64_t* group) {
	double* result = new double[6];
	for(int j = 0; j < N; j++) {
		result[group[j]]++;
	}
	return result;
}

double avg_op(int N, double* i) {
	double s = 0;
	for(int j = 0; j < N; j++) {
		s += i[j];
	}
	return s / N;
}

double* avg_op(int N, double* i, int64_t* group) {
	double* result = new double[6];
	int64_t* cnt = new int64_t[6];
	for(int j = 0; j < 6; j++) { result[j] = 0; cnt[j] = 0; }
	for(int j = 0; j < N; j++) {
		result[group[j]] += i[j];
		cnt[group[j]]++;
	}
	for(int j = 0; j < 6; j++) result[j] /= cnt[j];
	return result;
}

double sum_disc_price(int N, double* e, double* d) {
	double s = 0;
	for(int j = 0; j < N; j++) {
		s += e[j] * (1.0-d[j]);
	}
	return s;
}

double* sum_disc_price(int N, double* e, double* d, int64_t* group) {
	double* result = new double[6];
	for(int j = 0; j < 6; j++) result[j] = 0;
	for(int j = 0; j < N; j++) {
		result[group[j]] += e[j] * (1.0-d[j]);
	}
	return result;
}

double sum_charge(int N, double* e, double* d, double* t) {
	double s = 0;
	for(int j = 0; j < N; j++) {
		s += e[j] * (1.0-d[j]) * (1.0+t[j]);
	}
	return s;
}

double* sum_charge(int N, double* e, double* d, double* t, int64_t* group) {
	double* result = new double[6];
	for(int j = 0; j < 6; j++) result[j] = 0;
	for(int j = 0; j < N; j++) {
		result[group[j]] += e[j] * (1.0-d[j]) * (1.0+t[j]);
	}
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

	start_timing();

	int64_t* group = new int64_t[N];
	for(int i = 0; i < N; i++) {group[i] = (l_returnflag[i]+1) * (l_linestatus[i]+1) - 1;}

	double* sq = sum_op(N, l_quantity, group);
	double* sbp = sum_op(N, l_extendedprice, group);
	double* sdp = sum_disc_price(N, l_extendedprice, l_discount, group);
	double* sc = sum_charge(N, l_extendedprice, l_discount, l_tax, group);
	double* aq = avg_op(N, l_quantity, group);
	double* ap = avg_op(N, l_extendedprice, group);
	double* ad = avg_op(N, l_discount, group);
	double* count = count_op(N, group);

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
