# libdistance
Library for calculating distance between N-dimensional vectors. Supports all x86 SIMD instructions and ARM (aarch64) portable mode.

## Build and Test steps
### Step I. **Download google-test**
```bash
mkdir -p third-party
cd third-party
git clone https://github.com/google/googletest.git
cd ..
```
### Step II. **Build tests**
```bash
# It is important to build the library before the test
make clean all

# On x86_64 (Intel/AMD):
make test test_sse test_avx2 test_avx512
# On ARM (aarch64):
make test
```

### Step III. **Build Benchmark**
```bash
# On x86_64:
make bench bench_sse bench_avx2 bench_avx512
# On ARM (aarch64):
make bench
```

### Step IV. **Run Tests**
```bash
# Run the default (portable) test
./bin/test
# On x86_64, run SIMD tests:
./bin/test_sse
./bin/test_avx2
./bin/test_avx512
```

### Step V. **Run Benchmarks**
```bash
# Run the default (portable) benchmark
./bin/bench
# On x86_64, run SIMD benchmarks:
./bin/bench_sse
./bin/bench_avx2
./bin/bench_avx512
```

### Step VI. **Install Library and Header**
```bash
# Install default (portable) version
sudo make install
# On x86_64, install SIMD versions
sudo make install_sse
sudo make install_avx2
sudo make install_avx512
```
- Libraries are installed to `/usr/local/lib` as `libdistance.so` (and `libdistance_sse.so`, `libdistance_avx2.so`, `libdistance_avx512.so` on x86_64)
- Header is installed to `/usr/local/include/distance.h`

## API Reference

### Single Distance Functions
- `double euclidean_f64(const double* a, const double* b, size_t len);`
- `float  euclidean_f32(const float* a, const float* b, size_t len);`
- `double manhattan_f64(const double* a, const double* b, size_t len);`
- `float  manhattan_f32(const float* a, const float* b, size_t len);`

### Multi Distance Functions
- `double** multi_euclidean_f64(const double** a, const double** b, size_t len, size_t M, size_t N, size_t n_threads);`
- `float**  multi_euclidean_f32(const float** a, const float** b, size_t len, size_t M, size_t N, size_t n_threads);`
- `double** multi_manhattan_f64(const double** a, const double** b, size_t len, size_t M, size_t N, size_t n_threads);`
- `float**  multi_manhattan_f32(const float** a, const float** b, size_t len, size_t M, size_t N, size_t n_threads);`

## Features
- High-performance SIMD (SSE, AVX2, AVX512) and OpenMP parallelism (x86_64)
- Portable and optimized implementation for ARM (aarch64)
- Accurate and fast distance calculations for float and double
- Multi-vector (batch) distance computation
- GoogleTest-based unit tests for all functions
- Google Benchmark-based performance benchmarks
- Easy install and integration

## Requirements
- C++11 or later
- GCC or Clang with SIMD and OpenMP support
- [GoogleTest](https://github.com/google/googletest) and [Google Benchmark](https://github.com/google/benchmark) (for tests/benchmarks)

## Author
Siddhant Biradar