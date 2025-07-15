#include "distance.h"

double euclidean_f64(const double* vec_a, const double* vec_b, 
                     const size_t length){
#if defined (SSE)
    __m128d sum = _mm_setzero_pd();
    size_t i = 0;

    for(; i + 2 <= length; i += 2){
        __m128d va = _mm_loadu_pd(vec_a + i);
        __m128d vb = _mm_loadu_pd(vec_b + i);
        __m128d diff = _mm_sub_pd(va, vb);
        __m128d sq = _mm_mul_pd(diff, diff);
        sum = _mm_add_pd(sum, sq);
    }

    // Horizontal add for SSE
    __m128d high = _mm_unpackhi_pd(sum, sum);
    sum = _mm_add_sd(sum, high);
    double distance = _mm_cvtsd_f64(sum);

    // Handle remainder
    for(; i < length; i++) {
        double diff = vec_a[i] - vec_b[i];
        distance += diff * diff;
    }

    return sqrt(distance);

#elif defined (AVX2)
    __m256d sum = _mm256_setzero_pd();
    size_t i = 0;

    for(; i + 4 <= length; i += 4) {
        __m256d va = _mm256_loadu_pd(vec_a + i);
        __m256d vb = _mm256_loadu_pd(vec_b + i);
        __m256d diff = _mm256_sub_pd(va, vb);
        __m256d sq = _mm256_mul_pd(diff, diff);
        sum = _mm256_add_pd(sum, sq);
    }

    // Horizontal add for AVX2
    __m128d low = _mm256_castpd256_pd128(sum);
    __m128d high = _mm256_extractf128_pd(sum, 1);
    low = _mm_add_pd(low, high);
    high = _mm_unpackhi_pd(low, low);
    low = _mm_add_sd(low, high);
    double distance = _mm_cvtsd_f64(low);

    // Handle remainder
    for(; i < length; i++) {
        double diff = vec_a[i] - vec_b[i];
        distance += diff * diff;
    }

    return sqrt(distance);

#elif defined (AVX512)
    __m512d sum = _mm512_setzero_pd();
    size_t i = 0;

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
}

float euclidean_f32(const float* vec_a, const float* vec_b, 
                     const size_t length){
#if defined (SSE)
    __m128 sum = _mm_setzero_ps();
    size_t i = 0;

    for(; i + 4 <= length; i += 4){
        __m128 va = _mm_loadu_ps(vec_a + i);
        __m128 vb = _mm_loadu_ps(vec_b + i);
        __m128 diff = _mm_sub_ps(va, vb);
        __m128 sq = _mm_mul_ps(diff, diff);
        sum = _mm_add_ps(sum, sq);
    }

    // Horizontal add for SSE
    __m128 shuf = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(2, 3, 0, 1));
    sum = _mm_add_ps(sum, shuf);
    shuf = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(1, 0, 3, 2));
    sum = _mm_add_ps(sum, shuf);
    float distance = _mm_cvtss_f32(sum);

    // Handle remainder
    for(; i < length; i++) {
        float diff = vec_a[i] - vec_b[i];
        distance += diff * diff;
    }

    return sqrtf(distance);

#elif defined (AVX2)
    __m256 sum = _mm256_setzero_ps();
    size_t i = 0;

    for(; i + 8 <= length; i += 8) {
        __m256 va = _mm256_loadu_ps(vec_a + i);
        __m256 vb = _mm256_loadu_ps(vec_b + i);
        __m256 diff = _mm256_sub_ps(va, vb);
        __m256 sq = _mm256_mul_ps(diff, diff);
        sum = _mm256_add_ps(sum, sq);
    }

    // Horizontal add for AVX2
    __m128 low = _mm256_castps256_ps128(sum);
    __m128 high = _mm256_extractf128_ps(sum, 1);
    low = _mm_add_ps(low, high);
    __m128 shuf = _mm_shuffle_ps(low, low, _MM_SHUFFLE(2, 3, 0, 1));
    low = _mm_add_ps(low, shuf);
    shuf = _mm_shuffle_ps(low, low, _MM_SHUFFLE(1, 0, 3, 2));
    low = _mm_add_ps(low, shuf);
    float distance = _mm_cvtss_f32(low);

    // Handle remainder
    for(; i < length; i++) {
        float diff = vec_a[i] - vec_b[i];
        distance += diff * diff;
    }

    return sqrtf(distance);

#elif defined (AVX512)
    __m512 sum = _mm512_setzero_ps();
    size_t i = 0;

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
}

double manhattan_f64(const double* vec_a, const double* vec_b, 
                     const size_t length){
#if defined (SSE)
    __m128d sum = _mm_setzero_pd();
    size_t i = 0;

    for(; i + 2 <= length; i += 2){
        __m128d va = _mm_loadu_pd(vec_a + i);
        __m128d vb = _mm_loadu_pd(vec_b + i);
        __m128d diff = _mm_sub_pd(va, vb);
        __m128d abs_diff = _mm_andnot_pd(_mm_set1_pd(-0.0), diff);
        sum = _mm_add_pd(sum, abs_diff);
    }

    // Horizontal add for SSE
    __m128d high = _mm_unpackhi_pd(sum, sum);
    sum = _mm_add_sd(sum, high);
    double distance = _mm_cvtsd_f64(sum);

    // Handle remainder
    for(; i < length; i++) {
        double diff = vec_a[i] - vec_b[i];
        distance += fabs(diff);
    }

    return distance;

#elif defined (AVX2)
    __m256d sum = _mm256_setzero_pd();
    size_t i = 0;

    for(; i + 4 <= length; i += 4) {
        __m256d va = _mm256_loadu_pd(vec_a + i);
        __m256d vb = _mm256_loadu_pd(vec_b + i);
        __m256d diff = _mm256_sub_pd(va, vb);
        __m256d abs_diff = _mm256_andnot_pd(_mm256_set1_pd(-0.0), diff);
        sum = _mm256_add_pd(sum, abs_diff);
    }

    // Horizontal add for AVX2
    __m128d low = _mm256_castpd256_pd128(sum);
    __m128d high = _mm256_extractf128_pd(sum, 1);
    low = _mm_add_pd(low, high);
    high = _mm_unpackhi_pd(low, low);
    low = _mm_add_sd(low, high);
    double distance = _mm_cvtsd_f64(low);

    // Handle remainder
    for(; i < length; i++) {
        double diff = vec_a[i] - vec_b[i];
        distance += fabs(diff);
    }

    return distance;

#elif defined (AVX512)
    __m512d sum = _mm512_setzero_pd();
    size_t i = 0;

    for (; i + 8 <= length; i += 8) {
        __m512d va = _mm512_loadu_pd(vec_a + i);
        __m512d vb = _mm512_loadu_pd(vec_b + i);
        __m512d diff = _mm512_sub_pd(va, vb);
        __m512d abs_diff = _mm512_abs_pd(diff);
        sum = _mm512_add_pd(sum, abs_diff);
    }

    // Use reduce_add for horizontal sum
    double distance = _mm512_reduce_add_pd(sum);

    // Handle remainder
    for (; i < length; i++) {
        double diff = vec_a[i] - vec_b[i];
        distance += fabs(diff);
    }

    return distance;

#else
    double distance = 0.0;
    size_t i = 0;

    for(; i < length; i++) {
        double diff = vec_a[i] - vec_b[i];
        distance += fabs(diff);
    }

    return distance;
#endif
}

float manhattan_f32(const float* vec_a, const float* vec_b, 
                    const size_t length){
#if defined (SSE)
    __m128 sum = _mm_setzero_ps();
    size_t i = 0;

    for(; i + 4 <= length; i += 4){
        __m128 va = _mm_loadu_ps(vec_a + i);
        __m128 vb = _mm_loadu_ps(vec_b + i);
        __m128 diff = _mm_sub_ps(va, vb);
        __m128 abs_diff = _mm_andnot_ps(_mm_set1_ps(-0.0f), diff);
        sum = _mm_add_ps(sum, abs_diff);
    }

    // Horizontal add for SSE
    __m128 shuf = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(2, 3, 0, 1));
    sum = _mm_add_ps(sum, shuf);
    shuf = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(1, 0, 3, 2));
    sum = _mm_add_ps(sum, shuf);
    float distance = _mm_cvtss_f32(sum);

    // Handle remainder
    for(; i < length; i++) {
        float diff = vec_a[i] - vec_b[i];
        distance += fabsf(diff);
    }

    return distance;

#elif defined (AVX2)
    __m256 sum = _mm256_setzero_ps();
    size_t i = 0;

    for(; i + 8 <= length; i += 8) {
        __m256 va = _mm256_loadu_ps(vec_a + i);
        __m256 vb = _mm256_loadu_ps(vec_b + i);
        __m256 diff = _mm256_sub_ps(va, vb);
        __m256 abs_diff = _mm256_andnot_ps(_mm256_set1_ps(-0.0f), diff);
        sum = _mm256_add_ps(sum, abs_diff);
    }

    // Horizontal add for AVX2
    __m128 low = _mm256_castps256_ps128(sum);
    __m128 high = _mm256_extractf128_ps(sum, 1);
    low = _mm_add_ps(low, high);
    __m128 shuf = _mm_shuffle_ps(low, low, _MM_SHUFFLE(2, 3, 0, 1));
    low = _mm_add_ps(low, shuf);
    shuf = _mm_shuffle_ps(low, low, _MM_SHUFFLE(1, 0, 3, 2));
    low = _mm_add_ps(low, shuf);
    float distance = _mm_cvtss_f32(low);

    // Handle remainder
    for(; i < length; i++) {
        float diff = vec_a[i] - vec_b[i];
        distance += fabsf(diff);
    }

    return distance;

#elif defined (AVX512)
    __m512 sum = _mm512_setzero_ps();
    size_t i = 0;

    for (; i + 16 <= length; i += 16) {
        __m512 va = _mm512_loadu_ps(vec_a + i);
        __m512 vb = _mm512_loadu_ps(vec_b + i);
        __m512 diff = _mm512_sub_ps(va, vb);
        __m512 abs_diff = _mm512_abs_ps(diff);
        sum = _mm512_add_ps(sum, abs_diff);
    }

    // Use reduce_add for horizontal sum
    float distance = _mm512_reduce_add_ps(sum);

    // Handle remainder
    for (; i < length; i++) {
        float diff = vec_a[i] - vec_b[i];
        distance += fabsf(diff);
    }

    return distance;

#else
    float distance = 0.0;
    size_t i = 0;

    for(; i < length; i++) {
        float diff = vec_a[i] - vec_b[i];
        distance += fabsf(diff);
    }

    return distance;
#endif
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

double** multi_manhattan_f64(const double** vecs_a, const double** vecs_b, 
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
            distance_grid[i][j] = manhattan_f64(vecs_a[i], vecs_b[j], length);
        }
    }
    
    return distance_grid;
}

float** multi_manhattan_f32(const float** vecs_a, const float** vecs_b, 
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
    #pragma omp parallel for
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            distance_grid[i][j] = manhattan_f32(vecs_a[i], vecs_b[j], length);
        }
    }
    
    return distance_grid;
}