ldlib = -Wl,$(if find -name "*.so" $(LIBDIR),$(LIBDIR)/liblygemm.so,$(LIBDIR)/liblygemm.a) -lopenblas

test4x4: $(TESTDIR)/test4x4.cpp
	$(CXX) $(CXXFLAGS) -o $(TESTDIR)/test4x4 $< -L$(LIBDIR) $(ldlib)
	@echo "Running test"; 
	./$(TESTDIR)/test4x4; 

testclean:
	$(RM) $(TESTDIR)/test*[!.cpp]
