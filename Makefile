CC = gcc

CXX = g++
CFLAGS = -Wall -Wextra -fPIC -O3 -std=c11 -fopenmp
CXXFLAGS = -Wall -Wextra -O3 -std=c++17 -fopenmp
LDFLAGS = -lm

GTEST_SRC = third-party/googletest/googletest/src/gtest-all.cc
GTEST_MAIN_SRC = third-party/googletest/googletest/src/gtest_main.cc
GTEST_OBJS = $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o

BIN_DIR = bin


SRC = distance.c
OBJS = $(BIN_DIR)/libdistance.so $(BIN_DIR)/distance_sse.so $(BIN_DIR)/distance_avx2.so $(BIN_DIR)/distance_avx512.so

TEST_SRC = test.cc
TEST_BIN = $(BIN_DIR)/test
TEST_BIN_SSE = $(BIN_DIR)/test_sse
TEST_BIN_AVX2 = $(BIN_DIR)/test_avx2
TEST_BIN_AVX512 = $(BIN_DIR)/test_avx512

all: $(BIN_DIR) $(OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)


$(BIN_DIR)/libdistance.so: $(SRC)
	$(CC) $(CFLAGS) -shared $< -o $@ $(LDFLAGS)

$(BIN_DIR)/distance_sse.so: $(SRC)
	$(CC) $(CFLAGS) -DSSE -msse2 -shared $< -o $@ $(LDFLAGS)

$(BIN_DIR)/distance_avx2.so: $(SRC)
	$(CC) $(CFLAGS) -DAVX2 -mavx2 -mfma -shared $< -o $@ $(LDFLAGS)

$(BIN_DIR)/distance_avx512.so: $(SRC)
	$(CC) $(CFLAGS) -DAVX512 -mavx512f -shared $< -o $@ $(LDFLAGS)

gtest_main.o: $(GTEST_MAIN_SRC)

$(BIN_DIR)/gtest-all.o: $(GTEST_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest $< -c -o $@

$(BIN_DIR)/gtest_main.o: $(GTEST_MAIN_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest $< -c -o $@

test: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/libdistance.so
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/libdistance.so -o $(TEST_BIN)
	LD_LIBRARY_PATH=$(BIN_DIR) ./$(TEST_BIN)

test_sse: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/distance_sse.so
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/distance_sse.so -o $(TEST_BIN_SSE)
	LD_LIBRARY_PATH=$(BIN_DIR) LD_PRELOAD=$(BIN_DIR)/distance_sse.so $(TEST_BIN_SSE)

test_avx2: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/distance_avx2.so
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/distance_avx2.so -o $(TEST_BIN_AVX2)
	LD_LIBRARY_PATH=$(BIN_DIR) LD_PRELOAD=$(BIN_DIR)/distance_avx2.so $(TEST_BIN_AVX2)

test_avx512: $(GTEST_OBJS) $(TEST_SRC) $(BIN_DIR)/distance_avx512.so
	$(CXX) $(CXXFLAGS) -Ithird-party/googletest/googletest/include -Ithird-party/googletest/googletest -L$(BIN_DIR) -ldistance -lpthread -fopenmp -lm $(BIN_DIR)/gtest-all.o $(BIN_DIR)/gtest_main.o $(TEST_SRC) $(BIN_DIR)/distance_avx512.so -o $(TEST_BIN_AVX512)
	LD_LIBRARY_PATH=$(BIN_DIR) LD_PRELOAD=$(BIN_DIR)/distance_avx512.so $(TEST_BIN_AVX512)

clean:
	rm -f $(BIN_DIR)/*.so $(BIN_DIR)/*.o $(TEST_BIN)

.PHONY: all clean test
