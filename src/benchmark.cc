#include <benchmark/benchmark.h>

#include <cmath>
#include <vector>
#include <random>
#include <type_traits>

extern "C" {
#include "distance.h"
}

#define LENGTH 1024

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

static void bench_euclidean_f64(benchmark::State& state) {
    
    std::vector<double> A = random_vec<double>(LENGTH);
    std::vector<double> B = random_vec<double>(LENGTH);

    for (auto _ : state) {
        euclidean_f64(A.data(), B.data(), LENGTH);
    }
}
BENCHMARK(bench_euclidean_f64);

static void bench_euclidean_f32(benchmark::State& state) {
    
    std::vector<float> A = random_vec<float>(LENGTH);
    std::vector<float> B = random_vec<float>(LENGTH);

    for (auto _ : state) {
        euclidean_f32(A.data(), B.data(), LENGTH);
    }
}
BENCHMARK(bench_euclidean_f32);

static void bench_manhattan_f64(benchmark::State& state) {
    
    std::vector<double> A = random_vec<double>(LENGTH);
    std::vector<double> B = random_vec<double>(LENGTH);

    for (auto _ : state) {
        manhattan_f64(A.data(), B.data(), LENGTH);
    }
}
BENCHMARK(bench_manhattan_f64);

static void bench_manhattan_f32(benchmark::State& state) {
    
    std::vector<float> A = random_vec<float>(LENGTH);
    std::vector<float> B = random_vec<float>(LENGTH);

    for (auto _ : state) {
        manhattan_f32(A.data(), B.data(), LENGTH);
    }
}
BENCHMARK(bench_manhattan_f32);

static void bench_multi_euclidean_f64(benchmark::State& state) {
    constexpr size_t M = 32, N = 32;
    std::vector<std::vector<double>> vecs_a(M), vecs_b(N);
    for (size_t i = 0; i < M; ++i) vecs_a[i] = random_vec<double>(LENGTH);
    for (size_t j = 0; j < N; ++j) vecs_b[j] = random_vec<double>(LENGTH);
    std::vector<const double*> ptrs_a, ptrs_b;
    for (auto& v : vecs_a) ptrs_a.push_back(v.data());
    for (auto& v : vecs_b) ptrs_b.push_back(v.data());
    for (auto _ : state) {
        double** result = multi_euclidean_f64(ptrs_a.data(), ptrs_b.data(), LENGTH, M, N, 2);
        for (size_t i = 0; i < M; ++i) free(result[i]);
        free(result);
    }
}
BENCHMARK(bench_multi_euclidean_f64);

static void bench_multi_euclidean_f32(benchmark::State& state) {
    constexpr size_t M = 32, N = 32;
    std::vector<std::vector<float>> vecs_a(M), vecs_b(N);
    for (size_t i = 0; i < M; ++i) vecs_a[i] = random_vec<float>(LENGTH);
    for (size_t j = 0; j < N; ++j) vecs_b[j] = random_vec<float>(LENGTH);
    std::vector<const float*> ptrs_a, ptrs_b;
    for (auto& v : vecs_a) ptrs_a.push_back(v.data());
    for (auto& v : vecs_b) ptrs_b.push_back(v.data());
    for (auto _ : state) {
        float** result = multi_euclidean_f32(ptrs_a.data(), ptrs_b.data(), LENGTH, M, N, 2);
        for (size_t i = 0; i < M; ++i) free(result[i]);
        free(result);
    }
}
BENCHMARK(bench_multi_euclidean_f32);

static void bench_multi_manhattan_f64(benchmark::State& state) {
    constexpr size_t M = 32, N = 32;
    std::vector<std::vector<double>> vecs_a(M), vecs_b(N);
    for (size_t i = 0; i < M; ++i) vecs_a[i] = random_vec<double>(LENGTH);
    for (size_t j = 0; j < N; ++j) vecs_b[j] = random_vec<double>(LENGTH);
    std::vector<const double*> ptrs_a, ptrs_b;
    for (auto& v : vecs_a) ptrs_a.push_back(v.data());
    for (auto& v : vecs_b) ptrs_b.push_back(v.data());
    for (auto _ : state) {
        double** result = multi_manhattan_f64(ptrs_a.data(), ptrs_b.data(), LENGTH, M, N, 2);
        for (size_t i = 0; i < M; ++i) free(result[i]);
        free(result);
    }
}
BENCHMARK(bench_multi_manhattan_f64);

static void bench_multi_manhattan_f32(benchmark::State& state) {
    constexpr size_t M = 32, N = 32;
    std::vector<std::vector<float>> vecs_a(M), vecs_b(N);
    for (size_t i = 0; i < M; ++i) vecs_a[i] = random_vec<float>(LENGTH);
    for (size_t j = 0; j < N; ++j) vecs_b[j] = random_vec<float>(LENGTH);
    std::vector<const float*> ptrs_a, ptrs_b;
    for (auto& v : vecs_a) ptrs_a.push_back(v.data());
    for (auto& v : vecs_b) ptrs_b.push_back(v.data());
    for (auto _ : state) {
        float** result = multi_manhattan_f32(ptrs_a.data(), ptrs_b.data(), LENGTH, M, N, 2);
        for (size_t i = 0; i < M; ++i) free(result[i]);
        free(result);
    }
}
BENCHMARK(bench_multi_manhattan_f32);

static void bench_dot_product_f64(benchmark::State& state) {
    std::vector<double> vec_a = random_vec<double>(LENGTH);
    std::vector<double> vec_b = random_vec<double>(LENGTH);
    for (auto _ : state) {
        benchmark::DoNotOptimize(dot_product_f64(vec_a.data(), vec_b.data(), LENGTH));
    }
}
BENCHMARK(bench_dot_product_f64);

static void bench_dot_product_f32(benchmark::State& state) {
    std::vector<float> vec_a = random_vec<float>(LENGTH);
    std::vector<float> vec_b = random_vec<float>(LENGTH);
    for (auto _ : state) {
        benchmark::DoNotOptimize(dot_product_f32(vec_a.data(), vec_b.data(), LENGTH));
    }
}
BENCHMARK(bench_dot_product_f32);

static void bench_norm_f64(benchmark::State& state) {
    std::vector<double> vec = random_vec<double>(LENGTH);
    for (auto _ : state) {
        benchmark::DoNotOptimize(norm_f64(vec.data(), LENGTH));
    }
}
BENCHMARK(bench_norm_f64);

static void bench_norm_f32(benchmark::State& state) {
    std::vector<float> vec = random_vec<float>(LENGTH);
    for (auto _ : state) {
        benchmark::DoNotOptimize(norm_f32(vec.data(), LENGTH));
    }
}
BENCHMARK(bench_norm_f32);

static void bench_cosine_similarity_f64(benchmark::State& state) {
    std::vector<double> vec_a = random_vec<double>(LENGTH);
    std::vector<double> vec_b = random_vec<double>(LENGTH);
    for (auto _ : state) {
        benchmark::DoNotOptimize(cosine_similarity_f64(vec_a.data(), vec_b.data(), LENGTH));
    }
}
BENCHMARK(bench_cosine_similarity_f64);

static void bench_cosine_similarity_f32(benchmark::State& state) {
    std::vector<float> vec_a = random_vec<float>(LENGTH);
    std::vector<float> vec_b = random_vec<float>(LENGTH);
    for (auto _ : state) {
        benchmark::DoNotOptimize(cosine_similarity_f32(vec_a.data(), vec_b.data(), LENGTH));
    }
}
BENCHMARK(bench_cosine_similarity_f32);

BENCHMARK_MAIN();