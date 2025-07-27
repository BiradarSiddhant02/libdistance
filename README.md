# libdistance
High-performance library for calculating distances and mathematical operations between N-dimensional vectors. Supports x86 SIMD instructions (SSE3, AVX2, AVX512) with OpenMP parallelism and ARM (aarch64) NEON optimizations.

## Features
- **Distance Functions**: Euclidean, Manhattan, and **Minkowski** distance calculations with configurable p-norm
- **Mathematical Functions**: Dot product, vector norm (L2), and cosine similarity
- **High Performance**: SIMD-optimized implementations (SSE3, AVX2, AVX512) for x86_64
- **ARM Support**: NEON-optimized and portable fallback implementations for aarch64
- **Multi-threading**: OpenMP parallelization for batch operations
- **Comprehensive Testing**: GoogleTest unit tests with accuracy validation across all SIMD variants
- **Performance Benchmarking**: Google Benchmark integration for all architectures and distance metrics
- **Easy Integration**: Simple C API with installation support

## Directory Structure
```
libdistance/
├── include/distance.h     # Public API header
├── src/                   # Source files
│   ├── euclidean.c        # Euclidean distance implementations
│   ├── manhattan.c        # Manhattan distance implementations
│   ├── minkowski.c        # Minkowski distance implementations (generalized p-norm)
│   ├── cosine.c          # Dot product, norm, cosine similarity
│   ├── test.cc           # Comprehensive unit tests (16 test cases per architecture)
│   └── benchmark.cc      # Performance benchmarks (20+ benchmarks per architecture)
├── bin/                  # Build outputs
├── third-party/          # Dependencies (GoogleTest)
└── Makefile             # Build system

## Build and Test Steps

### Prerequisites
```bash
# Install dependencies (Ubuntu/Debian)
sudo apt update
sudo apt install build-essential gcc g++ libbenchmark-dev

# Clone GoogleTest (required for testing)
mkdir -p third-party
cd third-party
git clone https://github.com/google/googletest.git
cd ..
```

### Quick Start
```bash
# Build all libraries and run tests
make clean all
make test test_sse test_avx2  # Add test_avx512 if supported

# Run benchmarks
make bench bench_sse bench_avx2  # Add bench_avx512 if supported
```

### Detailed Build Process

#### Step I. **Clean and Build Libraries**
```bash
make clean    # Remove previous builds
make all      # Build all SIMD variants (default, SSE, AVX2, AVX512)
```

#### Step II. **Run Comprehensive Tests**
```bash
# On x86_64 (Intel/AMD):
make test         # Default portable implementation  
make test_sse     # SSE-optimized version
make test_avx2    # AVX2-optimized version  
make test_avx512  # AVX512-optimized version (if supported)

# On ARM (aarch64):
make test         # NEON-optimized version
```

#### Step III. **Run Performance Benchmarks**
```bash
# On x86_64:
make bench         # Default implementation
make bench_sse     # SSE benchmarks
make bench_avx2    # AVX2 benchmarks  
make bench_avx512  # AVX512 benchmarks (if supported)

# On ARM (aarch64):
make bench         # NEON benchmarks
```

#### Step IV. **Architecture Detection**
The Makefile automatically detects your system architecture:
```bash
make help         # View all available targets
make arch-info    # Display detected architecture and SIMD support
```

### Installation
```bash
# Install default (portable) version
sudo make install

# On x86_64, install specific SIMD versions
sudo make install_sse     # Install SSE version
sudo make install_avx2    # Install AVX2 version  
sudo make install_avx512  # Install AVX512 version (if supported)

# Uninstall
sudo make uninstall
```
- Libraries are installed to `/usr/local/lib/` 
- Header is installed to `/usr/local/include/distance.h`
- pkg-config files are provided for easy integration

## API Reference

### Distance Functions
```c
// Euclidean distance (L2 norm)
double euclidean_f64(const double* vec_a, const double* vec_b, size_t length);
float  euclidean_f32(const float* vec_a, const float* vec_b, size_t length);

// Manhattan distance (L1 norm)  
double manhattan_f64(const double* vec_a, const double* vec_b, size_t length);
float  manhattan_f32(const float* vec_a, const float* vec_b, size_t length);

// Minkowski distance (generalized Lp norm) - NEW!
double minkowski_f64(const double* vec_a, const double* vec_b, size_t length, double p);
float  minkowski_f32(const float* vec_a, const float* vec_b, size_t length, float p);
```

### Mathematical Functions
```c
// Dot product (inner product)
double dot_product_f64(const double* vec_a, const double* vec_b, size_t length);
float  dot_product_f32(const float* vec_a, const float* vec_b, size_t length);

// Vector norm (L2 norm, magnitude)
double norm_f64(const double* vec, size_t length);
float  norm_f32(const float* vec, size_t length);

// Cosine similarity  
double cosine_similarity_f64(const double* vec_a, const double* vec_b, size_t length);
float  cosine_similarity_f32(const float* vec_a, const float* vec_b, size_t length);
```

### Multi-Vector (Batch) Functions
```c
// Batch Euclidean distance computation
double** multi_euclidean_f64(const double** vecs_a, const double** vecs_b, 
                             size_t length, size_t M, size_t N, size_t n_threads);
float**  multi_euclidean_f32(const float** vecs_a, const float** vecs_b,
                             size_t length, size_t M, size_t N, size_t n_threads);

// Batch Manhattan distance computation
double** multi_manhattan_f64(const double** vecs_a, const double** vecs_b,
                             size_t length, size_t M, size_t N, size_t n_threads); 
float**  multi_manhattan_f32(const float** vecs_a, const float** vecs_b,
                             size_t length, size_t M, size_t N, size_t n_threads);
```

### Parameters
- `vec_a`, `vec_b`: Input vectors (aligned for optimal SIMD performance)
- `vec`: Input vector for single-vector operations
- `vecs_a`, `vecs_b`: Arrays of vector pointers for batch operations
- `length`: Vector dimension (number of elements)
- `p`: P-norm parameter for Minkowski distance (e.g., 1.0=Manhattan, 2.0=Euclidean, 0.5=fractional)
- `M`, `N`: Number of vectors in each batch (M × N distance matrix)
- `n_threads`: Number of OpenMP threads for parallelization

## Performance Characteristics

### SIMD Optimizations
- **Default**: Portable C implementation, works on all architectures
- **SSE3**: 128-bit SIMD with horizontal add instructions, processes 2 doubles or 4 floats per instruction
- **AVX2**: 256-bit SIMD with FMA support, processes 4 doubles or 8 floats per instruction  
- **AVX512**: 512-bit SIMD with reduce operations, processes 8 doubles or 16 floats per instruction
- **ARM NEON**: 128-bit SIMD optimizations for aarch64 with conditional ARMv8 features

### Multi-threading
- OpenMP parallelization for batch operations
- Automatic load balancing across available CPU cores
- Configurable thread count for optimal performance

## Example Usage

### Basic Distance Calculations
```c
#include <stdio.h>
#include <distance.h>

int main() {
    // Example vectors
    double vec_a[] = {1.0, 2.0, 3.0, 4.0};
    double vec_b[] = {5.0, 6.0, 7.0, 8.0};
    size_t length = 4;
    
    // Calculate distances
    double eucl_dist = euclidean_f64(vec_a, vec_b, length);
    double manh_dist = manhattan_f64(vec_a, vec_b, length);
    double mink_dist = minkowski_f64(vec_a, vec_b, length, 3.0);  // L3 norm
    
    // Calculate mathematical operations
    double dot_prod = dot_product_f64(vec_a, vec_b, length);
    double norm_a = norm_f64(vec_a, length);
    double cosine_sim = cosine_similarity_f64(vec_a, vec_b, length);
    
    printf("Euclidean distance: %f\n", eucl_dist);
    printf("Manhattan distance: %f\n", manh_dist);
    printf("Minkowski distance (p=3): %f\n", mink_dist);
    printf("Dot product: %f\n", dot_prod);
    printf("Norm of vec_a: %f\n", norm_a);
    printf("Cosine similarity: %f\n", cosine_sim);
    
    return 0;
}
```

### Batch Processing
```c
#include <stdlib.h>
#include <distance.h>

int main() {
    size_t M = 100, N = 50, length = 128, n_threads = 4;
    
    // Allocate vectors
    double** vecs_a = malloc(M * sizeof(double*));
    double** vecs_b = malloc(N * sizeof(double*));
    
    for (size_t i = 0; i < M; i++) {
        vecs_a[i] = malloc(length * sizeof(double));
        // Fill with data...
    }
    
    for (size_t j = 0; j < N; j++) {
        vecs_b[j] = malloc(length * sizeof(double));  
        // Fill with data...
    }
    
    // Compute M×N distance matrix in parallel
    double** distances = multi_euclidean_f64(
        (const double**)vecs_a, (const double**)vecs_b,
        length, M, N, n_threads
    );
    
    // Use distance matrix...
    printf("Distance[0][0]: %f\n", distances[0][0]);
    
    // Cleanup
    for (size_t i = 0; i < M; i++) {
        free(distances[i]);
        free(vecs_a[i]);
    }
    for (size_t j = 0; j < N; j++) {
        free(vecs_b[j]);
    }
    free(distances);
    free(vecs_a);
    free(vecs_b);
    
    return 0;
}
```

### Compilation
```bash
# Link with installed library
gcc -o example example.c -ldistance -lm

# Or use pkg-config  
gcc -o example example.c $(pkg-config --cflags --libs libdistance)
```

### Testing and Validation

### Unit Tests
The library includes comprehensive GoogleTest-based unit tests:
- **Accuracy Tests**: Validate correctness against known expected values
- **Edge Cases**: Handle zero vectors, single elements, large dimensions
- **Cross-Architecture**: Ensure consistent results across all SIMD implementations
- **Precision Tests**: Verify both float32 and float64 accuracy
- **Minkowski Tests**: Multiple p-values tested (0.5, 1.0, 1.5, 2.0, 3.0, 4.0) with appropriate tolerances

### Benchmarking
Google Benchmark integration provides:
- **Performance Measurement**: Throughput and latency analysis for all distance metrics
- **SIMD Comparison**: Performance gains across different instruction sets
- **Minkowski Benchmarks**: Performance testing across different p-values
- **Scaling Analysis**: Multi-threading efficiency evaluation
- **Architecture Profiling**: x86_64 vs ARM performance characteristics

### Benchmark Results Overview
The library includes 20+ benchmark cases per architecture:
- **Single-vector functions**: euclidean, manhattan, minkowski (5 p-values), dot_product, norm, cosine_similarity
- **Multi-vector functions**: batch euclidean and manhattan computations
- **Precision variants**: Both float32 and float64 for all functions
- **Minkowski p-value analysis**: Performance impact of different p-norm values (0.5, 1.0, 1.5, 2.0, 3.0)

## Advanced Features

### Minkowski Distance
The Minkowski distance is a generalization of many common distance metrics:
- **p = 0.5**: Fractional norm (computationally intensive)
- **p = 1.0**: Manhattan distance (L1 norm)
- **p = 1.5**: Custom intermediate norm
- **p = 2.0**: Euclidean distance (L2 norm)  
- **p = 3.0+**: Higher-order norms for specialized applications

All implementations use optimized horizontal sum intrinsics:
- **SSE3**: `_mm_hadd_ps/pd` for efficient vector reduction
- **AVX2**: `_mm256_hadd_ps/pd` with 256-bit operations
- **AVX512**: `_mm512_reduce_add_ps/pd` for maximum throughput
- **ARM NEON**: `vaddvq_f32/f64` where supported

### SIMD Implementation Details
- **Forced SSE3**: All SSE implementations use `hadd` instructions without conditional compilation
- **Horizontal Summation**: Optimized reduction operations for all architectures
- **Memory Alignment**: Proper vector alignment for maximum SIMD efficiency
- **Loop Unrolling**: Optimized inner loops for each SIMD variant

## Requirements
- **Compiler**: GCC 7+ or Clang 6+ with C11 and C++17 support  
- **SIMD Support**: SSE2+ (x86_64), NEON (ARM64) for optimized builds
- **OpenMP**: For multi-threading support (optional)
- **CMake**: 3.10+ (alternative build system)
- **Testing**: GoogleTest 1.10+ (for running tests)
- **Benchmarking**: Google Benchmark 1.5+ (for performance measurement)

## Contributing
1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Add tests for new functionality
4. Ensure all tests pass: `make test test_sse test_avx2`
5. Run benchmarks to verify performance: `make bench`
6. Submit a pull request

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Author
**Siddhant Biradar**  
Email: [your-email@example.com]  
GitHub: [@BiradarSiddhant02](https://github.com/BiradarSiddhant02)