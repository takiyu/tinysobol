#define TINY_SOBOL_IMPLEMENTATION
#include <cassert>
#include <set>
#include <sstream>

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

int main(int argc, char const* argv[]) {
#if 0
    // Single check
    CheckSobolNoOverlap(3, 1024);
#else
    // Check Sobol for each dimension and the number of samples
    for (uint32_t dim = 1; dim < 4; dim++) {
        std::cout << "Dim: " << dim << std::endl;
        uint32_t n_sample = 1;
        while (n_sample <= 1024) {
            std::cout << "  n_sample:" << n_sample << std::endl;
            const bool ret = CheckSobolNoOverlap<false>(dim, n_sample);
            if (ret) {
                std::cout << "    Passed" << std::endl;
            } else {
                std::cout << "    Has overlap" << std::endl;
            }
            n_sample <<= 1;
        }
    }
#endif

    return 0;
}
