#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <omp.h>

#include "distance.h"

double euclidean_f64(const double* vec_a, const double* vec_b, 
                     const size_t length){
#if defined (__x86_64__)
    #if defined (SSE)
        register __m128d sum = _mm_setzero_pd();
        register size_t i = 0;

        for(; i + 2 <= length; i += 2){
            __m128d va = _mm_loadu_pd(vec_a + i);
            __m128d vb = _mm_loadu_pd(vec_b + i);
            __m128d diff = _mm_sub_pd(va, vb);
            __m128d sq = _mm_mul_pd(diff, diff);
            sum = _mm_add_pd(sum, sq);
        }

        // Horizontal sum using SSE3 hadd instruction
        __m128d hadd_result = _mm_hadd_pd(sum, sum);
        double distance = _mm_cvtsd_f64(hadd_result);

        // Handle remainder
        for(; i < length; i++) {
            double diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrt(distance);

    #elif defined (AVX2)
        register __m256d sum = _mm256_setzero_pd();
        register size_t i = 0;

        for(; i + 4 <= length; i += 4) {
            __m256d va = _mm256_loadu_pd(vec_a + i);
            __m256d vb = _mm256_loadu_pd(vec_b + i);
            __m256d diff = _mm256_sub_pd(va, vb);
            __m256d sq = _mm256_mul_pd(diff, diff);
            sum = _mm256_add_pd(sum, sq);
        }

        // Horizontal sum using AVX hadd instructions
        __m256d hadd1 = _mm256_hadd_pd(sum, sum);
        __m128d sum_high = _mm256_extractf128_pd(hadd1, 1);
        __m128d sum_low = _mm256_castpd256_pd128(hadd1);
        __m128d final_sum = _mm_add_pd(sum_low, sum_high);
        double distance = _mm_cvtsd_f64(final_sum);

        // Handle remainder
        for(; i < length; i++) {
            double diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrt(distance);

    #elif defined (AVX512)
        register __m512d sum = _mm512_setzero_pd();
        register size_t i = 0;

        for (; i + 8 <= length; i += 8) {
            __m512d va = _mm512_loadu_pd(vec_a + i);
            __m512d vb = _mm512_loadu_pd(vec_b + i);
            __m512d diff = _mm512_sub_pd(va, vb);
            __m512d sq = _mm512_mul_pd(diff, diff);
            sum = _mm512_add_pd(sum, sq);
        }

        // Use reduce_add for horizontal sum
        double distance = _mm512_reduce_add_pd(sum);

        // Handle remainder
        for (; i < length; i++) {
            double diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrt(distance);

    #else
        double distance = 0.0;
        size_t i = 0;

        for(; i < length; i++) {
            double diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrt(distance);
    #endif
#elif defined (__aarch64__)
    register float64x2_t sum = vdupq_n_f64(0.0);  // Initialize vector to [0.0, 0.0]
    register size_t i = 0;

    // SIMD loop: process 2 doubles at a time
    for (; i + 2 <= length; i += 2) {
        float64x2_t va = vld1q_f64(vec_a + i);
        float64x2_t vb = vld1q_f64(vec_b + i);
        float64x2_t diff = vsubq_f64(va, vb);
        float64x2_t sq = vmulq_f64(diff, diff);
        sum = vaddq_f64(sum, sq);
    }

    // Horizontal sum - use ARMv8 vaddvq if available, otherwise manual
    #ifdef __ARM_ARCH_8A
        double distance = vaddvq_f64(sum);
    #else
        double distance = vgetq_lane_f64(sum, 0) + vgetq_lane_f64(sum, 1);
    #endif

    // Handle remaining elements (odd length)
    for (; i < length; ++i) {
        double diff = vec_a[i] - vec_b[i];
        distance += diff * diff;
    }

    return sqrt(distance);

#endif
    double distance1 = 0.0;
    size_t ii = 0;

    for(; ii < length; ii++) {
        double diff = vec_a[ii] - vec_b[ii];
        distance1 += diff * diff;
    }

    return sqrt(distance1);
}

float euclidean_f32(const float* vec_a, const float* vec_b, 
                     const size_t length){
#if defined (__x86_64__)
    #if defined (SSE)
        register __m128 sum = _mm_setzero_ps();
        register size_t i = 0;

        for(; i + 4 <= length; i += 4){
            __m128 va = _mm_loadu_ps(vec_a + i);
            __m128 vb = _mm_loadu_ps(vec_b + i);
            __m128 diff = _mm_sub_ps(va, vb);
            __m128 sq = _mm_mul_ps(diff, diff);
            sum = _mm_add_ps(sum, sq);
        }

        // Horizontal sum using SSE3 hadd instructions
        __m128 hadd1 = _mm_hadd_ps(sum, sum);
        __m128 hadd2 = _mm_hadd_ps(hadd1, hadd1);
        float distance = _mm_cvtss_f32(hadd2);

        // Handle remainder
        for(; i < length; i++) {
            float diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrtf(distance);

    #elif defined (AVX2)
        register __m256 sum = _mm256_setzero_ps();
        register size_t i = 0;

        for(; i + 8 <= length; i += 8) {
            __m256 va = _mm256_loadu_ps(vec_a + i);
            __m256 vb = _mm256_loadu_ps(vec_b + i);
            __m256 diff = _mm256_sub_ps(va, vb);
            __m256 sq = _mm256_mul_ps(diff, diff);
            sum = _mm256_add_ps(sum, sq);
        }

        // Horizontal sum using AVX hadd instructions
        __m256 hadd1 = _mm256_hadd_ps(sum, sum);
        __m256 hadd2 = _mm256_hadd_ps(hadd1, hadd1);
        __m128 sum_high = _mm256_extractf128_ps(hadd2, 1);
        __m128 sum_low = _mm256_castps256_ps128(hadd2);
        __m128 final_sum = _mm_add_ps(sum_low, sum_high);
        float distance = _mm_cvtss_f32(final_sum);

        // Handle remainder
        for(; i < length; i++) {
            float diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrtf(distance);

    #elif defined (AVX512)
        register __m512 sum = _mm512_setzero_ps();
        register size_t i = 0;

        for (; i + 16 <= length; i += 16) {
            __m512 va = _mm512_loadu_ps(vec_a + i);
            __m512 vb = _mm512_loadu_ps(vec_b + i);
            __m512 diff = _mm512_sub_ps(va, vb);
            __m512 sq = _mm512_mul_ps(diff, diff);
            sum = _mm512_add_ps(sum, sq);
        }

        // Use reduce_add for horizontal sum
        float distance = _mm512_reduce_add_ps(sum);

        // Handle remainder
        for (; i < length; i++) {
            float diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrtf(distance);

    #else
        float distance = 0.0;
        size_t i = 0;

        for(; i < length; i++) {
            float diff = vec_a[i] - vec_b[i];
            distance += diff * diff;
        }

        return sqrtf(distance);
    #endif
#elif defined (__aarch64__)
    register float32x4_t sum = vdupq_n_f32(0.0f);  // Initialize vector to [0.0, 0.0, 0.0, 0.0]
    register size_t i = 0;

    for (; i + 4 <= length; i += 4) {
        float32x4_t va = vld1q_f32(vec_a + i);
        float32x4_t vb = vld1q_f32(vec_b + i);
        float32x4_t diff = vsubq_f32(va, vb);
        float32x4_t sq = vmulq_f32(diff, diff);
        sum = vaddq_f32(sum, sq);
    }

    // Horizontal sum - use ARMv8 vaddvq if available, otherwise manual
    #ifdef __ARM_ARCH_8A
        float distance = vaddvq_f32(sum);
    #else
        float distance = vgetq_lane_f32(sum, 0) + vgetq_lane_f32(sum, 1)
                       + vgetq_lane_f32(sum, 2) + vgetq_lane_f32(sum, 3);
    #endif

    // Tail (scalar loop)
    for (; i < length; ++i) {
        float diff = vec_a[i] - vec_b[i];
        distance += diff * diff;
    }

    return sqrtf(distance);
#endif
    float distance1 = 0.0;
    size_t ii = 0;

    for(; ii < length; ii++) {
        float diff = vec_a[ii] - vec_b[ii];
        distance1 += diff * diff;
    }

    return sqrtf(distance1);
}

double** multi_euclidean_f64(const double** vecs_a, const double** vecs_b, 
                             const size_t length, const size_t M, const size_t N,
                             const size_t n_threads) {
    double** distance_grid = (double**)malloc(M * sizeof(double*));
    if (distance_grid == NULL) { 
        return NULL; 
    }
    
    for (size_t i = 0; i < M; i++) {
        distance_grid[i] = (double*)malloc(N * sizeof(double));
        if (distance_grid[i] == NULL) {
            for (size_t k = 0; k < i; k++) {
                free(distance_grid[k]);
            }
            free(distance_grid);
            return NULL;
        }
    }
    
    omp_set_num_threads(n_threads);
    #pragma omp parallel for
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            distance_grid[i][j] = euclidean_f64(vecs_a[i], vecs_b[j], length);
        }
    }
    
    return distance_grid;
}

float** multi_euclidean_f32(const float** vecs_a, const float** vecs_b, 
                            const size_t length, const size_t M, const size_t N,
                            const size_t n_threads) {
    float** distance_grid = (float**)malloc(M * sizeof(float*));
    if (distance_grid == NULL) { 
        return NULL; 
    }

    for (size_t i = 0; i < M; i++) {
        distance_grid[i] = (float*)malloc(N * sizeof(float));
        if (distance_grid[i] == NULL) {
            for (size_t k = 0; k < i; k++) {
                free(distance_grid[k]);
            }
            free(distance_grid);
            return NULL;
        }
    }

    omp_set_num_threads(n_threads);
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                distance_grid[i][j] = euclidean_f32(vecs_a[i], vecs_b[j], length);
            }
        }
    }
    
    return distance_grid;
}