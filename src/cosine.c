#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <omp.h>

#include "distance.h"

double dot_product_f64(const double* vec_a, const double* vec_b, 
                       const size_t length) {
#if defined (__x86_64__)
    #if defined (SSE)
        register __m128d sum = _mm_setzero_pd();
        register size_t i = 0;

        for (; i + 2 <= length; i += 2) {
            __m128d va = _mm_loadu_pd(vec_a + i);
            __m128d vb = _mm_loadu_pd(vec_b + i);
            __m128d prod = _mm_mul_pd(va, vb);
            sum = _mm_add_pd(sum, prod);
        }

        // Horizontal sum - use SSE3 hadd if available, otherwise manual
        #ifdef __SSE3__
            __m128d hadd_result = _mm_hadd_pd(sum, sum);
            double dot = _mm_cvtsd_f64(hadd_result);
        #else
            double buffer[2];
            _mm_store_pd(buffer, sum);
            double dot = buffer[0] + buffer[1];
        #endif

        for (; i < length; i++)
            dot += vec_a[i] * vec_b[i];

        return dot;

    #elif defined (AVX2)
        register __m256d sum = _mm256_setzero_pd();
        register size_t i = 0;

        for (; i + 4 <= length; i += 4) {
            __m256d va = _mm256_loadu_pd(vec_a + i);
            __m256d vb = _mm256_loadu_pd(vec_b + i);
            __m256d prod = _mm256_mul_pd(va, vb);
            sum = _mm256_add_pd(sum, prod);
        }

        // Horizontal sum using AVX hadd instructions
        __m256d hadd1 = _mm256_hadd_pd(sum, sum);
        __m128d sum_high = _mm256_extractf128_pd(hadd1, 1);
        __m128d sum_low = _mm256_castpd256_pd128(hadd1);
        __m128d final_sum = _mm_add_pd(sum_low, sum_high);
        double dot = _mm_cvtsd_f64(final_sum);

        for (; i < length; i++)
            dot += vec_a[i] * vec_b[i];

        return dot;

    #elif defined (AVX512)
        register __m512d sum = _mm512_setzero_pd();
        register size_t i = 0;
        
        for(; i + 8 <= length; i += 8) {
            __m512d va = _mm512_loadu_pd(vec_a + i);
            __m512d vb = _mm512_loadu_pd(vec_b + i);
            __m512d prod = _mm512_mul_pd(va, vb);
            sum = _mm512_add_pd(sum, prod);            
        }

        // Horizontal sum using AVX-512 reduce instruction
        double dot = _mm512_reduce_add_pd(sum);
        
        for (; i < length; i++)
            dot += vec_a[i] * vec_b[i];

        return dot;

    #else
        register double sum = 0.0;
        register size_t i = 0;

        for (; i < length; i++)
            sum += vec_a[i] * vec_b[i];
        
        return sum;

    #endif
#elif defined (__aarch64__)
    register float64x2_t sum = vdupq_n_f64(0.0);
    register size_t i = 0;

    for (; i + 2 <= length; i += 2) {
        float64x2_t va = vld1q_f64(vec_a + i);
        float64x2_t vb = vld1q_f64(vec_b + i);
        float64x2_t prod = vmulq_f64(va, vb);
        sum = vaddq_f64(sum, prod);
    }

    // Horizontal sum - use ARMv8 vaddvq if available, otherwise manual
    #ifdef __ARM_ARCH_8A
        double distance = vaddvq_f64(sum);
    #else
        double distance = vgetq_lane_f64(sum, 0) + vgetq_lane_f64(sum, 1);
    #endif

    for (; i < length; i++)
        distance += vec_a[i] * vec_b[i];
    
    return distance;

#else
    double prod1 = 0.0;

    for (size_t ii = 0; ii < length; ii++)
        prod1 += vec_a[ii] * vec_b[ii];

    return prod1;

#endif
}

float dot_product_f32(const float* vec_a, const float* vec_b, 
                       const size_t length) {
#if defined (__x86_64__)
    #if defined (SSE)
        register __m128 sum = _mm_setzero_ps();
        register size_t i = 0;

        for (; i + 4 <= length; i += 4) {
            __m128 va = _mm_loadu_ps(vec_a + i);
            __m128 vb = _mm_loadu_ps(vec_b + i);
            __m128 prod = _mm_mul_ps(va, vb);
            sum = _mm_add_ps(sum, prod);
        }

        // Horizontal sum - use SSE3 hadd if available, otherwise manual
        #ifdef __SSE3__
            __m128 hadd1 = _mm_hadd_ps(sum, sum);
            __m128 hadd2 = _mm_hadd_ps(hadd1, hadd1);
            float dot = _mm_cvtss_f32(hadd2);
        #else
            float buffer[4];
            _mm_store_ps(buffer, sum);
            float dot = buffer[0] + buffer[1] + buffer[2] + buffer[3];
        #endif

        for (; i < length; i++)
            dot += vec_a[i] * vec_b[i];

        return dot;

    #elif defined (AVX2)
        register __m256 sum = _mm256_setzero_ps();
        register size_t i = 0;

        for (; i + 8 <= length; i += 8) {
            __m256 va = _mm256_loadu_ps(vec_a + i);
            __m256 vb = _mm256_loadu_ps(vec_b + i);
            __m256 prod = _mm256_mul_ps(va, vb);
            sum = _mm256_add_ps(sum, prod);
        }

        // Horizontal sum using AVX hadd instructions
        __m256 hadd1 = _mm256_hadd_ps(sum, sum);
        __m256 hadd2 = _mm256_hadd_ps(hadd1, hadd1);
        __m128 sum_high = _mm256_extractf128_ps(hadd2, 1);
        __m128 sum_low = _mm256_castps256_ps128(hadd2);
        __m128 final_sum = _mm_add_ps(sum_low, sum_high);
        float dot = _mm_cvtss_f32(final_sum);

        for (; i < length; i++)
            dot += vec_a[i] * vec_b[i];

        return dot;

    #elif defined (AVX512)
        register __m512 sum = _mm512_setzero_ps();
        register size_t i = 0;
        
        for(; i + 16 <= length; i += 16) {
            __m512 va = _mm512_loadu_ps(vec_a + i);
            __m512 vb = _mm512_loadu_ps(vec_b + i);
            __m512 prod = _mm512_mul_ps(va, vb);
            sum = _mm512_add_ps(sum, prod);            
        }

        // Horizontal sum using AVX-512 reduce instruction
        float dot = _mm512_reduce_add_ps(sum);
        
        for (; i < length; i++)
            dot += vec_a[i] * vec_b[i];

        return dot;

    #else
        register float sum = 0.0f;
        register size_t i = 0;

        for (; i < length; i++)
            sum += vec_a[i] * vec_b[i];
        
        return sum;

    #endif
#elif defined (__aarch64__)
    register float32x4_t sum = vdupq_n_f32(0.0f);
    register size_t i = 0;

    for (; i + 4 <= length; i += 4) {
        float32x4_t va = vld1q_f32(vec_a + i);
        float32x4_t vb = vld1q_f32(vec_b + i);
        float32x4_t prod = vmulq_f32(va, vb);
        sum = vaddq_f32(sum, prod);
    }

    // Horizontal sum - use ARMv8 vaddvq if available, otherwise manual
    #ifdef __ARM_ARCH_8A
        float distance = vaddvq_f32(sum);
    #else
        float distance = vgetq_lane_f32(sum, 0) + vgetq_lane_f32(sum, 1)
                       + vgetq_lane_f32(sum, 2) + vgetq_lane_f32(sum, 3);
    #endif

    for (; i < length; i++)
        distance += vec_a[i] * vec_b[i];
    
    return distance;

#else
    float prod1 = 0.0f;

    for (size_t ii = 0; ii < length; ii++)
        prod1 += vec_a[ii] * vec_b[ii];

    return prod1;
    
#endif
}

double norm_f64(const double* vec, const size_t length) {
    return sqrt(dot_product_f64(vec, vec, length));
}

float norm_f32(const float* vec, const size_t length) {
    return sqrtf(dot_product_f32(vec, vec, length));
}

double cosine_similarity_f64(const double* vec_a, const double* vec_b, 
                             const size_t length) {
    return dot_product_f64(vec_a, vec_b, length) / 
                (norm_f64(vec_a, length) * norm_f64(vec_b, length));
}

float cosine_similarity_f32(const float* vec_a, const float* vec_b, 
                             const size_t length) {
    return dot_product_f32(vec_a, vec_b, length) / 
                (norm_f32(vec_a, length) * norm_f32(vec_b, length));
}