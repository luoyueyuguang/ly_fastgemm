# lygemm -- header-only small matrix multiplication
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -mavx2 -mfma -Wall -Wextra
INCLUDE  := -Iinclude
HEADERS  := $(wildcard include/*.hpp include/lygemm/*.hpp include/lygemm/detail/*.hpp)
BUILD    := build

# Comparison baselines for the benchmark, see bench/fetch_deps.sh
TP        := third_party
BENCH_INC := -Ibench -isystem $(TP)/root/usr/include -isystem $(TP)/root/usr/include/eigen3
BENCH_DEF := -DLYGEMM_BENCH_BLAS=1 -DLYGEMM_BENCH_XSMM=1
BENCH_LIB := -L$(TP)/lib -lxsmm -lopenblas -lpthread -ldl -lm
BENCH_RPATH := -Wl,-rpath,$(CURDIR)/$(TP)/root/usr/lib/x86_64-linux-gnu/openblas-pthread

.PHONY: test test-asan bench clean

test: $(BUILD)/test_correctness
	@$<

test-asan: $(BUILD)/test_correctness_asan
	@$<

$(BUILD)/test_correctness: tests/test_correctness.cpp $(HEADERS)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@

# the same suite under the sanitizers; the guarded buffers of the test catch
# out of bounds accesses that the sanitizers cannot see inside intrinsics
$(BUILD)/test_correctness_asan: tests/test_correctness.cpp $(HEADERS)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -fsanitize=address,undefined -fno-omit-frame-pointer $(INCLUDE) $< -o $@

bench: $(BUILD)/bench_gemm
	@$<

$(BUILD)/bench_gemm: bench/bench_gemm.cpp $(HEADERS) $(TP)/lib/libxsmm.a
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(BENCH_INC) $(BENCH_DEF) $< -o $@ $(BENCH_LIB) $(BENCH_RPATH)

$(TP)/lib/libxsmm.a:
	@echo "comparison baselines missing: run bench/fetch_deps.sh first" >&2
	@false

clean:
	$(RM) -r $(BUILD)
