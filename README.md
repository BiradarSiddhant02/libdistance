# libdistance
High-performance library for calculating distances and mathematical operations between N-dimensional vectors. Supports x86 SIMD instructions (SSE, AVX2, AVX512) with OpenMP parallelism and ARM (aarch64) NEON optimizations.

## Features
- **Distance Functions**: Euclidean and Manhattan distance calculations
- **Mathematical Functions**: Dot product, vector norm (L2), and cosine similarity
- **High Performance**: SIMD-optimized implementations (SSE, AVX2, AVX512) for x86_64
- **ARM Support**: NEON-optimized and portable fallback implementations for aarch64
- **Multi-threading**: OpenMP parallelization for batch operations
- **Comprehensive Testing**: GoogleTest unit tests with accuracy validation
- **Performance Benchmarking**: Google Benchmark integration for all architectures
- **Easy Integration**: Simple C API with installation support

## Directory Structure
```
libdistance/
├── include/distance.h     # Public API header
├── src/                   # Source files
│   ├── euclidean.c        # Euclidean distance implementations
│   ├── manhattan.c        # Manhattan distance implementations  
│   ├── cosine.c          # Dot product, norm, cosine similarity
│   ├── test.cc           # Comprehensive unit tests
│   └── benchmark.cc      # Performance benchmarks
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
- `M`, `N`: Number of vectors in each batch (M × N distance matrix)
- `n_threads`: Number of OpenMP threads for parallelization

## Performance Characteristics

### SIMD Optimizations
- **Default**: Portable C implementation, works on all architectures
- **SSE**: 128-bit SIMD, processes 2 doubles or 4 floats per instruction
- **AVX2**: 256-bit SIMD, processes 4 doubles or 8 floats per instruction  
- **AVX512**: 512-bit SIMD, processes 8 doubles or 16 floats per instruction
- **ARM NEON**: 128-bit SIMD optimizations for aarch64

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
    
    // Calculate mathematical operations
    double dot_prod = dot_product_f64(vec_a, vec_b, length);
    double norm_a = norm_f64(vec_a, length);
    double cosine_sim = cosine_similarity_f64(vec_a, vec_b, length);
    
    printf("Euclidean distance: %f\n", eucl_dist);
    printf("Manhattan distance: %f\n", manh_dist);
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

## Testing and Validation

### Unit Tests
The library includes comprehensive GoogleTest-based unit tests:
- **Accuracy Tests**: Validate correctness against known expected values
- **Edge Cases**: Handle zero vectors, single elements, large dimensions
- **Cross-Architecture**: Ensure consistent results across all SIMD implementations
- **Precision Tests**: Verify both float32 and float64 accuracy

### Benchmarking
Google Benchmark integration provides:
- **Performance Measurement**: Throughput and latency analysis
- **SIMD Comparison**: Performance gains across different instruction sets
- **Scaling Analysis**: Multi-threading efficiency evaluation
- **Architecture Profiling**: x86_64 vs ARM performance characteristics

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