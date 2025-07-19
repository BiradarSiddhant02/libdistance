
# Detect architecture
UNAME_ARCH := $(shell uname -m)
IS_X86_64 := $(filter x86_64,$(UNAME_ARCH))
IS_AARCH64 := $(filter aarch64,$(UNAME_ARCH))

CC = gcc
CXX = g++
CFLAGS = -Wall -Wextra -fPIC -O3 -std=c11 -fopenmp
CXXFLAGS = -Wall -Wextra -O3 -std=c++17 -fopenmp
LDFLAGS = -lm

GTEST_SRC = third-party/googletest/googletest/src/gtest-all.cc
GTEST_MAIN_SRC = third-party/googletest/googletest/src/gtest_main.cc
GTEST_OBJS = $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o

BIN_DIR = bin
SRC_DIR = src
INCLUDE_DIR = include


SRC = $(SRC_DIR)/euclidean.c $(SRC_DIR)/manhattan.c

ifeq ($(IS_X86_64),x86_64)
OBJS = $(BIN_DIR)/libdistance.so $(BIN_DIR)/libdistance_sse.so $(BIN_DIR)/libdistance_avx2.so $(BIN_DIR)/libdistance_avx512.so
else
OBJS = $(BIN_DIR)/libdistance.so
endif

TEST_SRC = $(SRC_DIR)/test.cc
TEST_BIN = $(BIN_DIR)/test
TEST_BIN_SSE = $(BIN_DIR)/test_sse
TEST_BIN_AVX2 = $(BIN_DIR)/test_avx2
TEST_BIN_AVX512 = $(BIN_DIR)/test_avx512


all: $(BIN_DIR) $(OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)



$(BIN_DIR)/libdistance.so: $(SRC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -shared $(SRC) -o $@ $(LDFLAGS)

ifeq ($(IS_X86_64),x86_64)
$(BIN_DIR)/libdistance_sse.so: $(SRC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -DSSE -msse2 -shared $(SRC) -o $@ $(LDFLAGS)

$(BIN_DIR)/libdistance_avx2.so: $(SRC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -DAVX2 -mavx2 -mfma -shared $(SRC) -o $@ $(LDFLAGS)

$(BIN_DIR)/libdistance_avx512.so: $(SRC)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -DAVX512 -mavx512f -shared $(SRC) -o $@ $(LDFLAGS)
endif

gtest_main.o: $(GTEST_MAIN_SRC)

$(BIN_DIR)/gtest-all.o: $(GTEST_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest $< -c -o $@

$(BIN_DIR)/gtest_main.o: $(GTEST_MAIN_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest $< -c -o $@


test: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/libdistance.so
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/libdistance.so -o $(TEST_BIN)
	LD_LIBRARY_PATH=$(BIN_DIR) ./$(TEST_BIN)

ifeq ($(IS_X86_64),x86_64)
test_sse: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/libdistance_sse.so
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/libdistance_sse.so -o $(TEST_BIN_SSE)
	LD_LIBRARY_PATH=$(BIN_DIR) LD_PRELOAD=$(BIN_DIR)/libdistance_sse.so $(TEST_BIN_SSE)

test_avx2: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/libdistance_avx2.so
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/libdistance_avx2.so -o $(TEST_BIN_AVX2)
	LD_LIBRARY_PATH=$(BIN_DIR) LD_PRELOAD=$(BIN_DIR)/libdistance_avx2.so $(TEST_BIN_AVX2)

test_avx512: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/libdistance_avx512.so
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/libdistance_avx512.so -o $(TEST_BIN_AVX512)
	LD_LIBRARY_PATH=$(BIN_DIR) LD_PRELOAD=$(BIN_DIR)/libdistance_avx512.so $(TEST_BIN_AVX512)
endif


install: install_default

install_default: $(BIN_DIR)/libdistance.so $(INCLUDE_DIR)/distance.h
	install -m 755 $(BIN_DIR)/libdistance.so /usr/local/lib/libdistance.so
	install -m 644 $(INCLUDE_DIR)/distance.h /usr/local/include/distance.h

ifeq ($(IS_X86_64),x86_64)
install_sse: $(BIN_DIR)/libdistance_sse.so $(INCLUDE_DIR)/distance.h
	install -m 755 $(BIN_DIR)/libdistance_sse.so /usr/local/lib/libdistance_sse.so
	install -m 644 $(INCLUDE_DIR)/distance.h /usr/local/include/distance.h

install_avx2: $(BIN_DIR)/libdistance_avx2.so $(INCLUDE_DIR)/distance.h
	install -m 755 $(BIN_DIR)/libdistance_avx2.so /usr/local/lib/libdistance_avx2.so
	install -m 644 $(INCLUDE_DIR)/distance.h /usr/local/include/distance.h

install_avx512: $(BIN_DIR)/libdistance_avx512.so $(INCLUDE_DIR)/distance.h
	install -m 755 $(BIN_DIR)/libdistance_avx512.so /usr/local/lib/libdistance_avx512.so
	install -m 644 $(INCLUDE_DIR)/distance.h /usr/local/include/distance.h
endif

# Usage:
#   make install           # installs default libdistance.so
#   make install_sse       # installs SSE version
#   make install_avx2      # installs AVX2 version
#   make install_avx512    # installs AVX512 version

BENCH_SRC = $(SRC_DIR)/benchmark.cc
BENCH_BIN = $(BIN_DIR)/bench
BENCH_BIN_SSE = $(BIN_DIR)/bench_sse
BENCH_BIN_AVX2 = $(BIN_DIR)/bench_avx2
BENCH_BIN_AVX512 = $(BIN_DIR)/bench_avx512


bench: $(BENCH_SRC) $(BIN_DIR)/libdistance.so
		$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -Ithird-party/benchmark/include -Lthird-party/benchmark/build/src -L$(BIN_DIR) \
			$(BENCH_SRC) -o $(BENCH_BIN) \
			-lbenchmark -ldistance -lpthread -fopenmp -lm
		LD_LIBRARY_PATH=$(BIN_DIR):third-party/benchmark/build/src ./$(BENCH_BIN)

ifeq ($(IS_X86_64),x86_64)
bench_sse: $(BENCH_SRC) $(BIN_DIR)/libdistance_sse.so
		$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -DSSE -msse2 -Ithird-party/benchmark/include -Lthird-party/benchmark/build/src -L$(BIN_DIR) \
			$(BENCH_SRC) -o $(BENCH_BIN_SSE) \
			-lbenchmark -ldistance -lpthread -fopenmp -lm
		LD_LIBRARY_PATH=$(BIN_DIR):third-party/benchmark/build/src LD_PRELOAD=$(BIN_DIR)/libdistance_sse.so ./$(BENCH_BIN_SSE)

bench_avx2: $(BENCH_SRC) $(BIN_DIR)/libdistance_avx2.so
		$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -DAVX2 -mavx2 -mfma -Ithird-party/benchmark/include -Lthird-party/benchmark/build/src -L$(BIN_DIR) \
			$(BENCH_SRC) -o $(BENCH_BIN_AVX2) \
			-lbenchmark -ldistance -lpthread -fopenmp -lm
		LD_LIBRARY_PATH=$(BIN_DIR):third-party/benchmark/build/src LD_PRELOAD=$(BIN_DIR)/libdistance_avx2.so ./$(BENCH_BIN_AVX2)

bench_avx512: $(BENCH_SRC) $(BIN_DIR)/libdistance_avx512.so
		$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -DAVX512 -mavx512f -Ithird-party/benchmark/include -Lthird-party/benchmark/build/src -L$(BIN_DIR) \
			$(BENCH_SRC) -o $(BENCH_BIN_AVX512) \
			-lbenchmark -ldistance -lpthread -fopenmp -lm
		LD_LIBRARY_PATH=$(BIN_DIR):third-party/benchmark/build/src LD_PRELOAD=$(BIN_DIR)/libdistance_avx512.so ./$(BENCH_BIN_AVX512)
endif


help:
	@echo "libdistance - High-performance vector distance library"
	@echo ""
	@echo "Architecture: $(UNAME_ARCH)"
ifeq ($(IS_X86_64),x86_64)
	@echo "SIMD support: SSE, AVX2, AVX512 available"
else ifeq ($(IS_AARCH64),aarch64)
	@echo "SIMD support: ARM NEON (portable mode)"
else
	@echo "SIMD support: Portable mode only"
endif
	@echo ""
	@echo "Available targets:"
	@echo "  all              - Build all libraries for this architecture"
	@echo "  test             - Build and run portable tests"
ifeq ($(IS_X86_64),x86_64)
	@echo "  test_sse         - Build and run SSE tests"
	@echo "  test_avx2        - Build and run AVX2 tests"
	@echo "  test_avx512      - Build and run AVX512 tests"
endif
	@echo "  bench            - Build and run portable benchmarks"
ifeq ($(IS_X86_64),x86_64)
	@echo "  bench_sse        - Build and run SSE benchmarks"
	@echo "  bench_avx2       - Build and run AVX2 benchmarks"
	@echo "  bench_avx512     - Build and run AVX512 benchmarks"
endif
	@echo "  install          - Install portable library and header"
ifeq ($(IS_X86_64),x86_64)
	@echo "  install_sse      - Install SSE library and header"
	@echo "  install_avx2     - Install AVX2 library and header"
	@echo "  install_avx512   - Install AVX512 library and header"
endif
	@echo "  clean            - Remove all build artifacts"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - GoogleTest: third-party/googletest/"
	@echo "  - Google Benchmark: third-party/benchmark/"

clean:
	rm -rf bin

.PHONY: all clean test help
