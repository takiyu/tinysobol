#define TINY_SOBOL_IMPLEMENTATION
#include <atomic>
#include <cassert>
#include <set>
#include <sstream>
#include <thread>

#include "tinysobol.h"

uint32_t Pow(uint32_t v, uint32_t n) {
    if (n == 1) {
        return v;
    } else {
        return Pow(v, n - 1) * v;
    }
}

void PrintSampledValue(uint32_t sample_idx, const std::vector<uint32_t>& sample,
                       bool is_overlapped) {
    std::cout << sample_idx << ": [ ";
    for (uint32_t i = 0; i < sample.size(); i++) {
        std::cout << sample[i] << " ";
    }
    std::cout << "]";
    if (is_overlapped) {
        std::cout << " Overlapped Value!";
    }
    std::cout << std::endl;
}

template <bool IsPrint = true>
bool CheckSobolNoOverlap(uint32_t dim, uint32_t n_sample,
                         uint32_t overlap_offset = 100) {
    // Create Sobol instance
    tinysobol::Sobol sobol(dim, n_sample);

    // Sampling loop
    std::set<std::vector<uint32_t>> sample_histroy;
    uint32_t overlapped_cnt = 0;
    for (uint32_t idx = 0; idx < Pow(n_sample, dim) + overlap_offset; idx++) {
        // Sample
        const std::vector<uint32_t>& sample = sobol.next();

        // Value check
        const bool is_overlapped = sample_histroy.count(sample);
        sample_histroy.emplace(sample);
        if (is_overlapped) {
            overlapped_cnt++;
        }

        // Value print
        if (IsPrint) {
            PrintSampledValue(idx, sample, is_overlapped);
        }
    }

    // Check overlapped count
    return overlapped_cnt == overlap_offset;
}

void PrintDimNSample(uint32_t dim, uint32_t n_sample, const std::string& msg) {
    std::cout << "Dim:" << dim << ", n_sample:" << n_sample << "  >> " << msg
              << std::endl;
}

int main(int argc, char const* argv[]) {
#if 0
    // Single check
    CheckSobolNoOverlap(3, 1024);
#else
    // Check all combinations
    const uint32_t N_THREADS = std::thread::hardware_concurrency();
    for (uint32_t dim = 1; dim < 4; dim++) {
        std::atomic<uint32_t> n_sample = {1};
        std::vector<std::thread> threads;
        for (uint32_t t_idx = 0; t_idx < N_THREADS; t_idx++) {
            threads.emplace_back([&]() {
                while (n_sample++ <= 1024) {
                    const bool ret = CheckSobolNoOverlap<false>(dim, n_sample);
                    if (ret) {
                        PrintDimNSample(dim, n_sample, "Passed");
                    } else {
                        PrintDimNSample(dim, n_sample, "Has invalid overlap");
                    }
                }
            });
        }
        for (auto&& thread : threads) {
            thread.join();
        }
    }
#endif

    return 0;
}
