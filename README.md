# libdistance
Library for calculating distance between N-dimensional vectors. Supports all x86 SIMD instructions

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

# Build tests
make test test_sse test_avx2 test_avx512
```
