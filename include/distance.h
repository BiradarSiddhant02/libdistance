#ifndef DISTANCE_H
#define DISTANCE_H

#if defined(__x86_64__)
    #if defined(SSE) || defined(AVX2) || defined(AVX512)
        #include <immintrin.h>
    #endif
#elif defined(__aarch64__)
    #include <arm_neon.h>
#else
    #warning "SIMD not enabled or unsupported architecture."
#endif

#if defined (__cplusplus)
extern "C"{
#endif
/** 
 * @brief Function to calculate euclidean distance between two 
 * vectors at 64 bit precision.
 *
 * @param vec_a (const double*) vector A
 * @param vec_b (const double*) vector B
 * @param length (const size_t) Length of vectors
 *
 * @return distance between Vector A and Vector B
 */
double euclidean_f64(const double*, const double*, const size_t);

/** 
 * @brief Function to calculate euclidean distance between two
 * vectors at 32 bit precision.
 *
 * @param vec_a (const float*) Vector A
 * @param vec_b (const float*) Vector B
 * @param length (const size_t) Length of vectors
 *
 * @return distance between Vector A and Vector B
 */
float euclidean_f32(const float*, const float*, const size_t);

/** 
 * @brief Function to calculate manhattan distance between two
 * vectors at 64 bit precision.
 *
 * @param vec_a (const double*) Vector A
 * @param vec_b (const double*) Vector B
 * @param length (const size_t) Length of vectors
 *
 * @return distance between Vector A and Vector B
 */
double manhattan_f64(const double*, const double*, const size_t);

/** 
 * @brief Function to calculate manhattan distance between two
 * vectors at 32 bit precision.
 *
 * @param vec_a (const float*) Vector A
 * @param vec_b (const float*) Vector B
 * @param length (const size_t) Length of vectors
 *
 * @return distance between Vector A and Vector B
 */
float manhattan_f32(const float*, const float*, const size_t);

/**
 * @brief Function to calculate dot product between two
 * vectors in 64 bit precision.
 * 
 * @param vec_a (const double*) Vector A
 * @param vec_b (const double*) Vector B
 * @param length (const size_t) Length of vectors
 * 
 * @return Dot product of two vectors
 */
double dot_product_f64(const double*, const double*, const size_t);

/**
 * @brief Function to calculate dot product between two
 * vectors in 32 bit precision.
 * 
 * @param vec_a (const float*) Vector A
 * @param vec_b (const float*) Vector B
 * @param length (const size_t) Length of vectors
 * 
 * @return Dot product of two vectors
 */
float dot_product_f32(const float*, const float*, const size_t);

/**
 * @brief Function to calculate norm of a vector at
 * 64 bit precision
 * 
 * @param vec (const double*) Vector
 * @param length (const size_t) length of vector
 * 
 * @return norm of the vector
 */
double norm_f64(const double*, const size_t);

/**
 * @brief Function to calculate norm of a vector at
 * 32 bit precision
 * 
 * @param vec (const float*) Vector
 * @param length (const size_t) length of vector
 * 
 * @return norm of the vector
 */
float norm_f32(const float*, const size_t);

/**
 * @brief Function to calculate cosine similarity at
 * 64 bit precision
 * 
 * @param vec_a (const double*) Vector A
 * @param vec_b (const double*) Vector B
 * @param length (const size_t) length of both vectors
 * 
 * @return Value between -1 and 1
 */
double cosine_similarity_f64(const double*, const double*, const size_t);

/**
 * @brief Function to calculate cosine similarity at
 * 32 bit precision
 * 
 * @param vec_a (const float*) Vector A
 * @param vec_b (const float*) Vector B
 * @param length (const size_t) length of both vectors
 * 
 * @return Value between -1 and 1
 */
float cosine_similarity_f32(const float*, const float*, const size_t);

/**
 * @brief Function to calculate euclidean distance between M and N vectors
 * with 64 bit precision
 * 
 * @param vecs_a (const double**) M Vectors
 * @param vecs_b (const double**) N Vectors
 * @param length Length of each vector
 * @param M Number of vectors in vecs_a
 * @param N Number of vectors in vecs_b
 * @param n_threads Number of threads
 * 
 * @return 2D grid of distance between M and N vectors
 */
double** multi_euclidean_f64(const double**, const double**, const size_t,
                             const size_t, const size_t, const size_t);

/**
 * @brief Function to calculate euclidean distance between M and N vectors
 * with 32 bit precision
 * 
 * @param vecs_a (const double**) M Vectors
 * @param vecs_b (const double**) N Vectors
 * @param length Length of each vector
 * @param M Number of vectors in vecs_a
 * @param N Number of vectors in vecs_b
 * @param n_threads Number of threads
 * 
 * @return 2D grid of distance between M and N vectors
 */
float** multi_euclidean_f32(const float**, const float**, const size_t,
                             const size_t, const size_t, const size_t);

/**
 * @brief Function to calculate manhattan distance between M and N vectors
 * with 64 bit precision
 * 
 * @param vecs_a (const double**) M Vectors
 * @param vecs_b (const double**) N Vectors
 * @param length Length of each vector
 * @param M Number of vectors in vecs_a
 * @param N Number of vectors in vecs_b
 * @param n_threads Number of threads
 * 
 * @return 2D grid of distance between M and N vectors
 */
double** multi_manhattan_f64(const double**, const double**, const size_t,
                             const size_t, const size_t, const size_t);

/**
 * @brief Function to calculate manhattan distance between M and N vectors
 * with 32 bit precision
 * 
 * @param vecs_a (const double**) M Vectors
 * @param vecs_b (const double**) N Vectors
 * @param length Length of each vector
 * @param M Number of vectors in vecs_a
 * @param N Number of vectors in vecs_b
 * @param n_threads Number of threads
 * 
 * @return 2D grid of distance between M and N vectors
 */
float** multi_manhattan_f32(const float**, const float**, const size_t,
                            const size_t, const size_t, const size_t);
#if defined (__cplusplus)
}
#endif

#endif // DISTANCE_H