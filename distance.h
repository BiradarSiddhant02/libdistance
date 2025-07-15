#include <immintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <math.h>
#include <omp.h>

/** 
 * @brief Function to calculate euclidean distance between two 
 * vectors at 64 bit precision.
 *
 * @param vec_a (const double*) vector A
 * @param vec_b (const double*) vector B
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
 *
 * @return distance between Vector A and Vector B
 */
float manhattan_f32(const float*, const float*, const size_t);

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