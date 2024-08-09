#!/usr/bin/bash

CXXFLAGS += -I include

VPATH += src include lib

$(LIBDIR)/liblygemm.a(lygemm.o): lygemm.o

$(LIBDIR)/liblygemm.so: lygemm.o 
	$(CXX) -shared -o $(@) $(^)

libclean:
	$(RM) $(LIBDIR)/*.o $(LIBDIR)/*.a $(LIBDIR)/*.so
