CPP=c++
CPPFLAGS=-std=c++17 -Wall -Wextra -Werror -O2

.PHONY: clean

out-main: cmain.o dft.o fft.o ssa.o
	env $(CPP) $^ $(CPPFLAGS) -o $@

cmain.o: cmain.cpp
	env $(CPP) $< $(CPPFLAGS) -c -o $@

dft.o: dft.cpp dft.hpp
	env $(CPP) $< $(CPPFLAGS) -c -o $@

fft.o: fft.cpp fft.hpp
	env $(CPP) $< $(CPPFLAGS) -c -o $@

ssa.o: ssa.cpp ssa.hpp
	env $(CPP) $< $(CPPFLAGS) -c -o $@

clean:
	env rm -rf out-main *.o
