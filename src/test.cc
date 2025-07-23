#include <gtest/gtest.h>

#include <cmath>
#include <vector>
#include <random>
#include <type_traits>

template <typename T>
class RandomGenerator {
    static_assert(std::is_floating_point<T>::value, "Template type must be float or double");

public:
    RandomGenerator(T min, T max)
        : dist(min, max), rng(std::random_device{}()) {}

    T operator()() {
        return dist(rng);
    }

private:
    std::conditional_t<std::is_same<T, float>::value,
                       std::uniform_real_distribution<float>,
                       std::uniform_real_distribution<double>> dist;
    std::mt19937 rng;
};

static RandomGenerator<float> rng_f32(-1.0f, 1.0f);
static RandomGenerator<double> rng_f64(-1.0d, 1.0d);

template <typename precT>
std::vector<precT> random_vec(std::size_t length) {
    std::vector<precT> vec(length);
    for (std::size_t i = 0; i < length; ++i)
        vec[i] = std::is_same<precT, float>::value ? rng_f32() : rng_f64();
    return vec;
}

extern "C" {
#include "distance.h"
}

template <typename precT>
precT expected_euclidean(const std::vector<precT>& A, const std::vector<precT>& B) {
    precT distance = static_cast<precT>(0.0);
    for (size_t i = 0; i < A.size(); i++) {
        precT diff = A[i] - B[i];
        distance += diff * diff;
    }
    return static_cast<precT>(std::sqrt(distance));
}

template <typename precT>
precT expected_manhattan(const std::vector<precT>& A, const std::vector<precT>& B) {
    precT distance = static_cast<precT>(0.0);
    for (size_t i = 0; i < A.size(); i++) {
        precT diff = A[i] - B[i];
        distance += std::abs(diff);
    }
    return distance;
}

template <typename precT>
precT expected_dot_product(const std::vector<precT>& A, const std::vector<precT>& B) {
    precT dot_product = static_cast<precT>(0.0);
    for (size_t i = 0; i < A.size(); i++) {
        dot_product += A[i] * B[i];
    }
    return dot_product;
}

template <typename precT>
precT expected_norm(const std::vector<precT>& A) {
    precT norm = static_cast<precT>(0.0);
    for (size_t i = 0; i < A.size(); i++) {
        norm += A[i] * A[i];
    }
    return static_cast<precT>(std::sqrt(norm));
}

template <typename precT>
precT expected_cosine_similarity(const std::vector<precT>& A, const std::vector<precT>& B) {
    precT dot_product = expected_dot_product(A, B);
    precT norm_a = expected_norm(A);
    precT norm_b = expected_norm(B);
    if (norm_a == static_cast<precT>(0.0) || norm_b == static_cast<precT>(0.0)) {
        return static_cast<precT>(0.0);
    }
    return dot_product / (norm_a * norm_b);
}

template <typename precT>
void test_euclidean_accuracy() {
    static const size_t length = 128;
    std::vector<precT> A = random_vec<precT>(length);
    std::vector<precT> B = random_vec<precT>(length);
    precT expected = expected_euclidean(A, B);
    precT actual;
    if constexpr (std::is_same<precT, float>::value) {
        actual = euclidean_f32(A.data(), B.data(), length);
    } else {
        actual = euclidean_f64(A.data(), B.data(), length);
    }
    ASSERT_NEAR(expected, actual, static_cast<precT>(1e-4));
}

template <typename precT>
void test_manhattan_accuracy() {
    static const size_t length = 128;
    std::vector<precT> A = random_vec<precT>(length);
    std::vector<precT> B = random_vec<precT>(length);
    precT expected = expected_manhattan(A, B);
    precT actual;
    if constexpr (std::is_same<precT, float>::value) {
        actual = manhattan_f32(A.data(), B.data(), length);
    } else {
        actual = manhattan_f64(A.data(), B.data(), length);
    }
    ASSERT_NEAR(expected, actual, static_cast<precT>(1e-4));
}

template <typename precT>
void test_dot_product_accuracy() {
    static const size_t length = 128;
    std::vector<precT> A = random_vec<precT>(length);
    std::vector<precT> B = random_vec<precT>(length);
    precT expected = expected_dot_product(A, B);
    precT actual;
    if constexpr (std::is_same<precT, float>::value) {
        actual = dot_product_f32(A.data(), B.data(), length);
    } else {
        actual = dot_product_f64(A.data(), B.data(), length);
    }
    ASSERT_NEAR(expected, actual, static_cast<precT>(1e-4));
}

template <typename precT>
void test_norm_accuracy() {
    static const size_t length = 128;
    std::vector<precT> A = random_vec<precT>(length);
    precT expected = expected_norm(A);
    precT actual;
    if constexpr (std::is_same<precT, float>::value) {
        actual = norm_f32(A.data(), length);
    } else {
        actual = norm_f64(A.data(), length);
    }
    ASSERT_NEAR(expected, actual, static_cast<precT>(1e-4));
}

template <typename precT>
void test_cosine_similarity_accuracy() {
    static const size_t length = 128;
    std::vector<precT> A = random_vec<precT>(length);
    std::vector<precT> B = random_vec<precT>(length);
    precT expected = expected_cosine_similarity(A, B);
    precT actual;
    if constexpr (std::is_same<precT, float>::value) {
        actual = cosine_similarity_f32(A.data(), B.data(), length);
    } else {
        actual = cosine_similarity_f64(A.data(), B.data(), length);
    }
    ASSERT_NEAR(expected, actual, static_cast<precT>(1e-4));
}

template <typename precT>
void test_multi_euclidean_accuracy() {
    static const size_t length = 64;
    static const size_t M = 8, N = 8;
    std::vector<std::vector<precT>> vecs_a(M), vecs_b(N);
    for (size_t i = 0; i < M; ++i) vecs_a[i] = random_vec<precT>(length);
    for (size_t j = 0; j < N; ++j) vecs_b[j] = random_vec<precT>(length);
    std::vector<const precT*> ptrs_a, ptrs_b;
    for (auto& v : vecs_a) ptrs_a.push_back(v.data());
    for (auto& v : vecs_b) ptrs_b.push_back(v.data());
    precT** result;
    if constexpr (std::is_same<precT, float>::value) {
        result = multi_euclidean_f32(ptrs_a.data(), ptrs_b.data(), length, M, N, 2);
    } else {
        result = multi_euclidean_f64(ptrs_a.data(), ptrs_b.data(), length, M, N, 2);
    }
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            precT expected = expected_euclidean(vecs_a[i], vecs_b[j]);
            ASSERT_NEAR(expected, result[i][j], static_cast<precT>(1e-4));
        }
    }
    for (size_t i = 0; i < M; ++i) free(result[i]);
    free(result);
}

template <typename precT>
void test_multi_manhattan_accuracy() {
    static const size_t length = 64;
    static const size_t M = 8, N = 8;
    std::vector<std::vector<precT>> vecs_a(M), vecs_b(N);
    for (size_t i = 0; i < M; ++i) vecs_a[i] = random_vec<precT>(length);
    for (size_t j = 0; j < N; ++j) vecs_b[j] = random_vec<precT>(length);
    std::vector<const precT*> ptrs_a, ptrs_b;
    for (auto& v : vecs_a) ptrs_a.push_back(v.data());
    for (auto& v : vecs_b) ptrs_b.push_back(v.data());
    precT** result;
    if constexpr (std::is_same<precT, float>::value) {
        result = multi_manhattan_f32(ptrs_a.data(), ptrs_b.data(), length, M, N, 2);
    } else {
        result = multi_manhattan_f64(ptrs_a.data(), ptrs_b.data(), length, M, N, 2);
    }
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            precT expected = expected_manhattan(vecs_a[i], vecs_b[j]);
            ASSERT_NEAR(expected, result[i][j], static_cast<precT>(1e-4));
        }
    }
    for (size_t i = 0; i < M; ++i) free(result[i]);
    free(result);
}

// --- TESTS FOR ALL COMPILATION CASES ---

#if defined(SSE)
TEST(TestDistanceSSE, EuclideanF32) { test_euclidean_accuracy<float>(); }
TEST(TestDistanceSSE, EuclideanF64) { test_euclidean_accuracy<double>(); }
TEST(TestDistanceSSE, ManhattanF32) { test_manhattan_accuracy<float>(); }
TEST(TestDistanceSSE, ManhattanF64) { test_manhattan_accuracy<double>(); }
TEST(TestDistanceSSE, DotProductF32) { test_dot_product_accuracy<float>(); }
TEST(TestDistanceSSE, DotProductF64) { test_dot_product_accuracy<double>(); }
TEST(TestDistanceSSE, NormF32) { test_norm_accuracy<float>(); }
TEST(TestDistanceSSE, NormF64) { test_norm_accuracy<double>(); }
TEST(TestDistanceSSE, CosineSimilarityF32) { test_cosine_similarity_accuracy<float>(); }
TEST(TestDistanceSSE, CosineSimilarityF64) { test_cosine_similarity_accuracy<double>(); }
TEST(TestDistanceSSE, MultiEuclideanF32) { test_multi_euclidean_accuracy<float>(); }
TEST(TestDistanceSSE, MultiEuclideanF64) { test_multi_euclidean_accuracy<double>(); }
TEST(TestDistanceSSE, MultiManhattanF32) { test_multi_manhattan_accuracy<float>(); }
TEST(TestDistanceSSE, MultiManhattanF64) { test_multi_manhattan_accuracy<double>(); }
#elif defined(AVX2)
TEST(TestDistanceAVX2, EuclideanF32) { test_euclidean_accuracy<float>(); }
TEST(TestDistanceAVX2, EuclideanF64) { test_euclidean_accuracy<double>(); }
TEST(TestDistanceAVX2, ManhattanF32) { test_manhattan_accuracy<float>(); }
TEST(TestDistanceAVX2, ManhattanF64) { test_manhattan_accuracy<double>(); }
TEST(TestDistanceAVX2, DotProductF32) { test_dot_product_accuracy<float>(); }
TEST(TestDistanceAVX2, DotProductF64) { test_dot_product_accuracy<double>(); }
TEST(TestDistanceAVX2, NormF32) { test_norm_accuracy<float>(); }
TEST(TestDistanceAVX2, NormF64) { test_norm_accuracy<double>(); }
TEST(TestDistanceAVX2, CosineSimilarityF32) { test_cosine_similarity_accuracy<float>(); }
TEST(TestDistanceAVX2, CosineSimilarityF64) { test_cosine_similarity_accuracy<double>(); }
TEST(TestDistanceAVX2, MultiEuclideanF32) { test_multi_euclidean_accuracy<float>(); }
TEST(TestDistanceAVX2, MultiEuclideanF64) { test_multi_euclidean_accuracy<double>(); }
TEST(TestDistanceAVX2, MultiManhattanF32) { test_multi_manhattan_accuracy<float>(); }
TEST(TestDistanceAVX2, MultiManhattanF64) { test_multi_manhattan_accuracy<double>(); }
#elif defined(AVX512)
TEST(TestDistanceAVX512, EuclideanF32) { test_euclidean_accuracy<float>(); }
TEST(TestDistanceAVX512, EuclideanF64) { test_euclidean_accuracy<double>(); }
TEST(TestDistanceAVX512, ManhattanF32) { test_manhattan_accuracy<float>(); }
TEST(TestDistanceAVX512, ManhattanF64) { test_manhattan_accuracy<double>(); }
TEST(TestDistanceAVX512, DotProductF32) { test_dot_product_accuracy<float>(); }
TEST(TestDistanceAVX512, DotProductF64) { test_dot_product_accuracy<double>(); }
TEST(TestDistanceAVX512, NormF32) { test_norm_accuracy<float>(); }
TEST(TestDistanceAVX512, NormF64) { test_norm_accuracy<double>(); }
TEST(TestDistanceAVX512, CosineSimilarityF32) { test_cosine_similarity_accuracy<float>(); }
TEST(TestDistanceAVX512, CosineSimilarityF64) { test_cosine_similarity_accuracy<double>(); }
TEST(TestDistanceAVX512, MultiEuclideanF32) { test_multi_euclidean_accuracy<float>(); }
TEST(TestDistanceAVX512, MultiEuclideanF64) { test_multi_euclidean_accuracy<double>(); }
TEST(TestDistanceAVX512, MultiManhattanF32) { test_multi_manhattan_accuracy<float>(); }
TEST(TestDistanceAVX512, MultiManhattanF64) { test_multi_manhattan_accuracy<double>(); }
#else
TEST(TestDistanceDefault, EuclideanF32) { test_euclidean_accuracy<float>(); }
TEST(TestDistanceDefault, EuclideanF64) { test_euclidean_accuracy<double>(); }
TEST(TestDistanceDefault, ManhattanF32) { test_manhattan_accuracy<float>(); }
TEST(TestDistanceDefault, ManhattanF64) { test_manhattan_accuracy<double>(); }
TEST(TestDistanceDefault, DotProductF32) { test_dot_product_accuracy<float>(); }
TEST(TestDistanceDefault, DotProductF64) { test_dot_product_accuracy<double>(); }
TEST(TestDistanceDefault, NormF32) { test_norm_accuracy<float>(); }
TEST(TestDistanceDefault, NormF64) { test_norm_accuracy<double>(); }
TEST(TestDistanceDefault, CosineSimilarityF32) { test_cosine_similarity_accuracy<float>(); }
TEST(TestDistanceDefault, CosineSimilarityF64) { test_cosine_similarity_accuracy<double>(); }
TEST(TestDistanceDefault, MultiEuclideanF32) { test_multi_euclidean_accuracy<float>(); }
TEST(TestDistanceDefault, MultiEuclideanF64) { test_multi_euclidean_accuracy<double>(); }
TEST(TestDistanceDefault, MultiManhattanF32) { test_multi_manhattan_accuracy<float>(); }
TEST(TestDistanceDefault, MultiManhattanF64) { test_multi_manhattan_accuracy<double>(); }
#endif
