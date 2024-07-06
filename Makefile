#!/usr/bin/bash
all:
	$(CXX) -c -fPIC src/* -Iinclude -mavx2 -mfma -O3 
	if [ ! -d lib ]; then mkdir lib; fi
	$(CXX) -shared -o lib/liblygemm.so lygemm.o 
	mv lygemm.o lib/
	rm -f src/lygemm.hpp.gch

test:
	g++ -o tests/test4x4 tests/test4x4.cpp -Llib -llygemm -mavx2 -mfma -O3 -Iinclude -lopenblas
	LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH ./tests/test4x4


clean:
	rm -f lib/liblygemm.so lib/lygemm.o 
	if [ -f tests/test4x4 ]; then rm tests/test4x4; fi
