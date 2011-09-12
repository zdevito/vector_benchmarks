GCCFLAGS=-O3 -msse4 -DNDEBUG -funroll-loops -DROUNDS=20 -DLENGTH=1048576
GCCFLAGS += $(FLAGS)
GCC=g++

TESTS = fused fused_sse strip_mined strip_mined_sse long_to_short_interpreter

all:	GCCFLAGS += -DBLOCK=64
all:	$(TESTS)

test:	GCCFLAGS += $(FLAGS)
test:	$(TESTS)

asm:	GCCFLAGS += -S -DBLOCK=64
asm:	$(TESTS)

fused:	fused.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
fused_sse:	fused_sse.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
fused_avx:	fused_avx.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
strip_mined:	strip_mined.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
strip_mined_sse:	strip_mined_sse.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
strip_mined_avx:	strip_mined_avx.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
long_to_short_interpreter:	long_to_short_interpreter.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
clean:
	rm $(TESTS)
