#define TINY_SOBOL_IMPLEMENTATION
#include <atomic>
#include <cassert>
#include <set>
#include <sstream>
#include <thread>

#include "tinysobol.h"

uint64_t Pow(uint64_t v, uint64_t n) {
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
    for (uint64_t idx = 0; idx < Pow(n_sample, dim) + overlap_offset; idx++) {
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
    bool ret = CheckSobolNoOverlap(3, 1024 * 2);
    std::cout << ret << std::endl;
#else
    // Check all combinations
    const uint32_t N_THREADS = std::thread::hardware_concurrency();
    const uint32_t N_DIM_MAX = 10;
    const uint32_t N_SAMPLE_SHIFT_MAX = 12;
    std::atomic<uint32_t> idx_pool = {1};
    std::vector<std::thread> threads;
    for (uint32_t t_idx = 0; t_idx < N_THREADS; t_idx++) {
        threads.emplace_back([&]() {
            uint32_t idx = 0;
            while ((idx = idx_pool++) <= N_DIM_MAX * N_SAMPLE_SHIFT_MAX) {
                uint32_t dim = (idx / N_SAMPLE_SHIFT_MAX) + 1;
                uint32_t n_sample_shift = idx % N_SAMPLE_SHIFT_MAX;
                uint32_t n_sample = 1 << n_sample_shift;
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
#endif

    return 0;
}
