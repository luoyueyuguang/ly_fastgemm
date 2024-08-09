CXX ?= icpx
CXX ?= g++
CXX ?= clang++

CXXFLAGS += -O3 -mfma -mavx2 

LIBDIR=lib
INCDIR=include
TESTDIR=tests

.PHONY: all static shared clean test
.INTERMEDIATE: lygemm.o

all:static shared

static:$(LIBDIR)/liblygemm.a(lygemm.o)

shared:$(LIBDIR)/liblygemm.so

test:test4x4

clean:libclean inclean testclean
	
include lib/Makefile include/Makefile tests/Makefile 
