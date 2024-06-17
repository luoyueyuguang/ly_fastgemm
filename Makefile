#!/usr/bin/bash
CXX = clang++
all:
	$(CXX) -c -fPIC src/lygemm.cpp src/lygemm.hpp -mavx2 -mfma -O3 
	$(CXX) -shared -o lib/liblygemm.so lygemm.o -lopenblas
	mv lygemm.o lib/
	rm -f src/lygemm.hpp.gch

clean:
	rm -f lib/liblygemm.so lib/lygemm.o
