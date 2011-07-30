#include <stdint.h>
__inline__ uint64_t rdtsc() {
	uint32_t low, high;
	__asm__ __volatile__ (
		"xorl %%eax,%%eax \n    cpuid"
		::: "%rax", "%rbx", "%rcx", "%rdx" );
	__asm__ __volatile__ (
						  "rdtsc" : "=a" (low), "=d" (high));
	return (uint64_t)high << 32 | low;
}





#include <sys/time.h>

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

