#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <omp.h>

#include "distance.h"

double minkowski_f64(const double* vec_a, const double* vec_b, 
                     const size_t length, const double p) {
#if defined (__x86_64__)
    #if defined (SSE)
        register __m128d sum = _mm_setzero_pd();
        register size_t i = 0;

        for (; i + 2 <= length; i += 2) {
            __m128d va = _mm_loadu_pd(vec_a + i);
            __m128d vb = _mm_loadu_pd(vec_b + i);
            __m128d diff = _mm_sub_pd(va, vb);
            __m128d abs_diff = _mm_andnot_pd(_mm_set1_pd(-0.0), diff);
            
            // Calculate abs_diff^p for each element
            double temp[2];
            _mm_store_pd(temp, abs_diff);
            temp[0] = pow(temp[0], p);
            temp[1] = pow(temp[1], p);
            __m128d powered = _mm_load_pd(temp);
            
            sum = _mm_add_pd(sum, powered);
        }

        // Horizontal sum using SSE3 hadd instruction
        __m128d hadd_result = _mm_hadd_pd(sum, sum);
        double distance = _mm_cvtsd_f64(hadd_result);

        // Handle remainder
        for (; i < length; i++) {
            double diff = fabs(vec_a[i] - vec_b[i]);
            distance += pow(diff, p);
        }

        return pow(distance, 1.0 / p);

    #elif defined (AVX2)
        register __m256d sum = _mm256_setzero_pd();
        register size_t i = 0;

        for (; i + 4 <= length; i += 4) {
            __m256d va = _mm256_loadu_pd(vec_a + i);
            __m256d vb = _mm256_loadu_pd(vec_b + i);
            __m256d diff = _mm256_sub_pd(va, vb);
            __m256d abs_diff = _mm256_andnot_pd(_mm256_set1_pd(-0.0), diff);
            
            // Calculate abs_diff^p for each element
            double temp[4];
            _mm256_store_pd(temp, abs_diff);
            temp[0] = pow(temp[0], p);
            temp[1] = pow(temp[1], p);
            temp[2] = pow(temp[2], p);
            temp[3] = pow(temp[3], p);
            __m256d powered = _mm256_load_pd(temp);
            
            sum = _mm256_add_pd(sum, powered);
        }

        // Horizontal sum using AVX hadd instructions
        __m256d hadd1 = _mm256_hadd_pd(sum, sum);
        __m128d sum_high = _mm256_extractf128_pd(hadd1, 1);
        __m128d sum_low = _mm256_castpd256_pd128(hadd1);
        __m128d final_sum = _mm_add_pd(sum_low, sum_high);
        double distance = _mm_cvtsd_f64(final_sum);

        // Handle remainder
        for (; i < length; i++) {
            double diff = fabs(vec_a[i] - vec_b[i]);
            distance += pow(diff, p);
        }

        return pow(distance, 1.0 / p);

    #elif defined (AVX512)
        register __m512d sum = _mm512_setzero_pd();
        register size_t i = 0;
        
        for (; i + 8 <= length; i += 8) {
            __m512d va = _mm512_loadu_pd(vec_a + i);
            __m512d vb = _mm512_loadu_pd(vec_b + i);
            __m512d diff = _mm512_sub_pd(va, vb);
            __m512d abs_diff = _mm512_abs_pd(diff);
            
            // Calculate abs_diff^p for each element
            double temp[8];
            _mm512_store_pd(temp, abs_diff);
            for (int j = 0; j < 8; j++) {
                temp[j] = pow(temp[j], p);
            }
            __m512d powered = _mm512_load_pd(temp);
            
            sum = _mm512_add_pd(sum, powered);
        }

        // Horizontal sum using AVX-512 reduce instruction
        double distance = _mm512_reduce_add_pd(sum);
        
        // Handle remainder
        for (; i < length; i++) {
            double diff = fabs(vec_a[i] - vec_b[i]);
            distance += pow(diff, p);
        }

        return pow(distance, 1.0 / p);

    #else
        register double sum = 0.0;
        register size_t i = 0;

        for (; i < length; i++) {
            double diff = fabs(vec_a[i] - vec_b[i]);
            sum += pow(diff, p);
        }
        
        return pow(sum, 1.0 / p);

    #endif
#elif defined (__aarch64__)
    register float64x2_t sum = vdupq_n_f64(0.0);
    register size_t i = 0;

    for (; i + 2 <= length; i += 2) {
        float64x2_t va = vld1q_f64(vec_a + i);
        float64x2_t vb = vld1q_f64(vec_b + i);
        float64x2_t diff = vsubq_f64(va, vb);
        float64x2_t abs_diff = vabsq_f64(diff);
        
        // Calculate abs_diff^p for each element
        double temp[2];
        vst1q_f64(temp, abs_diff);
        temp[0] = pow(temp[0], p);
        temp[1] = pow(temp[1], p);
        float64x2_t powered = vld1q_f64(temp);
        
        sum = vaddq_f64(sum, powered);
    }

    // Horizontal sum - use ARMv8 vaddvq if available, otherwise manual
    #ifdef __ARM_ARCH_8A
        double distance = vaddvq_f64(sum);
    #else
        double distance = vgetq_lane_f64(sum, 0) + vgetq_lane_f64(sum, 1);
    #endif

    for (; i < length; i++) {
        double diff = fabs(vec_a[i] - vec_b[i]);
        distance += pow(diff, p);
    }
    
    return pow(distance, 1.0 / p);

#else
    double sum = 0.0;

    for (size_t ii = 0; ii < length; ii++) {
        double diff = fabs(vec_a[ii] - vec_b[ii]);
        sum += pow(diff, p);
    }

    return pow(sum, 1.0 / p);

#endif
}

float minkowski_f32(const float* vec_a, const float* vec_b, 
                    const size_t length, const float p) {
#if defined (__x86_64__)
    #if defined (SSE)
        register __m128 sum = _mm_setzero_ps();
        register size_t i = 0;

        for (; i + 4 <= length; i += 4) {
            __m128 va = _mm_loadu_ps(vec_a + i);
            __m128 vb = _mm_loadu_ps(vec_b + i);
            __m128 diff = _mm_sub_ps(va, vb);
            __m128 abs_diff = _mm_andnot_ps(_mm_set1_ps(-0.0f), diff);
            
            // Calculate abs_diff^p for each element
            float temp[4];
            _mm_store_ps(temp, abs_diff);
            temp[0] = powf(temp[0], p);
            temp[1] = powf(temp[1], p);
            temp[2] = powf(temp[2], p);
            temp[3] = powf(temp[3], p);
            __m128 powered = _mm_load_ps(temp);
            
            sum = _mm_add_ps(sum, powered);
        }

        // Horizontal sum using SSE3 hadd instructions
        __m128 hadd1 = _mm_hadd_ps(sum, sum);
        __m128 hadd2 = _mm_hadd_ps(hadd1, hadd1);
        float distance = _mm_cvtss_f32(hadd2);

        // Handle remainder
        for (; i < length; i++) {
            float diff = fabsf(vec_a[i] - vec_b[i]);
            distance += powf(diff, p);
        }

        return powf(distance, 1.0f / p);

    #elif defined (AVX2)
        register __m256 sum = _mm256_setzero_ps();
        register size_t i = 0;

        for (; i + 8 <= length; i += 8) {
            __m256 va = _mm256_loadu_ps(vec_a + i);
            __m256 vb = _mm256_loadu_ps(vec_b + i);
            __m256 diff = _mm256_sub_ps(va, vb);
            __m256 abs_diff = _mm256_andnot_ps(_mm256_set1_ps(-0.0f), diff);
            
            // Calculate abs_diff^p for each element
            float temp[8];
            _mm256_store_ps(temp, abs_diff);
            for (int j = 0; j < 8; j++) {
                temp[j] = powf(temp[j], p);
            }
            __m256 powered = _mm256_load_ps(temp);
            
            sum = _mm256_add_ps(sum, powered);
        }

        // Horizontal sum using AVX hadd instructions
        __m256 hadd1 = _mm256_hadd_ps(sum, sum);
        __m256 hadd2 = _mm256_hadd_ps(hadd1, hadd1);
        __m128 sum_high = _mm256_extractf128_ps(hadd2, 1);
        __m128 sum_low = _mm256_castps256_ps128(hadd2);
        __m128 final_sum = _mm_add_ps(sum_low, sum_high);
        float distance = _mm_cvtss_f32(final_sum);

        // Handle remainder
        for (; i < length; i++) {
            float diff = fabsf(vec_a[i] - vec_b[i]);
            distance += powf(diff, p);
        }

        return powf(distance, 1.0f / p);

    #elif defined (AVX512)
        register __m512 sum = _mm512_setzero_ps();
        register size_t i = 0;
        
        for (; i + 16 <= length; i += 16) {
            __m512 va = _mm512_loadu_ps(vec_a + i);
            __m512 vb = _mm512_loadu_ps(vec_b + i);
            __m512 diff = _mm512_sub_ps(va, vb);
            __m512 abs_diff = _mm512_abs_ps(diff);
            
            // Calculate abs_diff^p for each element
            float temp[16];
            _mm512_store_ps(temp, abs_diff);
            for (int j = 0; j < 16; j++) {
                temp[j] = powf(temp[j], p);
            }
            __m512 powered = _mm512_load_ps(temp);
            
            sum = _mm512_add_ps(sum, powered);
        }

        // Horizontal sum using AVX-512 reduce instruction
        float distance = _mm512_reduce_add_ps(sum);
        
        // Handle remainder
        for (; i < length; i++) {
            float diff = fabsf(vec_a[i] - vec_b[i]);
            distance += powf(diff, p);
        }

        return powf(distance, 1.0f / p);

    #else
        register float sum = 0.0f;
        register size_t i = 0;

        for (; i < length; i++) {
            float diff = fabsf(vec_a[i] - vec_b[i]);
            sum += powf(diff, p);
        }
        
        return powf(sum, 1.0f / p);

    #endif
#elif defined (__aarch64__)
    register float32x4_t sum = vdupq_n_f32(0.0f);
    register size_t i = 0;

    for (; i + 4 <= length; i += 4) {
        float32x4_t va = vld1q_f32(vec_a + i);
        float32x4_t vb = vld1q_f32(vec_b + i);
        float32x4_t diff = vsubq_f32(va, vb);
        float32x4_t abs_diff = vabsq_f32(diff);
        
        // Calculate abs_diff^p for each element
        float temp[4];
        vst1q_f32(temp, abs_diff);
        temp[0] = powf(temp[0], p);
        temp[1] = powf(temp[1], p);
        temp[2] = powf(temp[2], p);
        temp[3] = powf(temp[3], p);
        float32x4_t powered = vld1q_f32(temp);
        
        sum = vaddq_f32(sum, powered);
    }

    // Horizontal sum - use ARMv8 vaddvq if available, otherwise manual
    #ifdef __ARM_ARCH_8A
        float distance = vaddvq_f32(sum);
    #else
        float distance = vgetq_lane_f32(sum, 0) + vgetq_lane_f32(sum, 1)
                       + vgetq_lane_f32(sum, 2) + vgetq_lane_f32(sum, 3);
    #endif

    for (; i < length; i++) {
        float diff = fabsf(vec_a[i] - vec_b[i]);
        distance += powf(diff, p);
    }
    
    return powf(distance, 1.0f / p);

#else
    float sum = 0.0f;

    for (size_t ii = 0; ii < length; ii++) {
        float diff = fabsf(vec_a[ii] - vec_b[ii]);
        sum += powf(diff, p);
    }

    return powf(sum, 1.0f / p);
    
#endif
}