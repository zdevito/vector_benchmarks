GCCFLAGS=-O3 -DNDEBUG -funroll-loops
GCC=g++

TESTS = interpreter strip_mined compiled long_interpreter

all:	$(TESTS)

interpreter:	interpreter.cpp
	$(GCC) $(GCCFLAGS) $< -o $@

strip_mined:	strip_mined.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
compiled:	compiled.cpp
	$(GCC) $(GCCFLAGS) $< -o $@
	
long_interpreter:	long_interpreter.cpp
	$(GCC) $(GCCFLAGS) $< -o $@

	
clean:
	rm $(TESTS)