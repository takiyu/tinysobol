/*
  License:
    MIT license.

  Notation:
    This code is based on `https://github.com/naught101/sobol_seq`.

*/

#ifndef TINY_SOBOL_H_190529
#define TINY_SOBOL_H_190529

#include <array>
#include <cassert>
#include <cstdint>
#include <exception>
#include <iostream>
#include <vector>

namespace tinysobol {

// -----------------------------------------------------------------------------
// ----------------------------------- Sobol -----------------------------------
// -----------------------------------------------------------------------------
static constexpr uint64_t DIM_MAX = 40;
static constexpr uint64_t LOG_MAX = 30;
using VGrid = std::array<std::array<uint64_t, LOG_MAX>, DIM_MAX>;  // [DIM][LOG]

class Sobol {
public:
    Sobol(uint64_t dim, uint64_t n_sample, uint64_t seed = 0);
    Sobol(uint64_t dim, const std::vector<uint64_t>& n_samples,
          uint64_t seed = 0);

    const std::vector<uint64_t>& next();
    const std::vector<uint64_t>& getLastSample() const;

private:
    const uint64_t DIM;
    const VGrid V_GRID;
    const uint64_t RECIPD_DENOM;
    const std::vector<uint64_t> N_SAMPLES_2BASE;
    const std::vector<uint64_t> N_SAMPLES;

    uint64_t m_seed;
    std::vector<uint64_t> m_last_q;
    std::vector<uint64_t> m_last_sample;

    void nextImpl();
    bool isValidSample();
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// --------------------------- Begin of Definitions ----------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
#ifdef TINY_SOBOL_IMPLEMENTATION

namespace {

uint64_t GetHighestBitPos(uint64_t n) {
    // Returns the position of the high 1 bit base 2 in an integer.
    uint64_t bit = 0;
    while (n > 0) {
        bit++;
        n >>= 1;
    }
    return bit;
}

uint64_t GetLowestBitPos(uint64_t n) {
    // Returns the position of the low 0 bit base 2 in an integer.
    uint64_t bit = 1;
    while (n != (n & (~uint64_t(1)))) {
        bit++;
        n >>= 1;
    }
    return bit;
}

constexpr uint64_t V_TMPL[DIM_MAX][LOG_MAX] = {
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 3, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 7, 13, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 5, 11, 27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 5, 1, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 7, 3, 29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 7, 7, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 1, 9, 23, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 3, 5, 19, 33, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 3, 13, 11, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 7, 13, 25, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 5, 11, 7, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 1, 3, 13, 39, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 1, 15, 17, 63, 13, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 5, 5, 1, 27, 33, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 3, 3, 25, 17, 115, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0,   0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 3, 15, 29, 15, 41, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 1, 7, 3, 23, 79, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 7, 9, 31, 29, 17, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 5, 13, 11, 3, 29, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 1, 9, 5, 21, 119, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,  0,   0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 3, 1, 23, 13, 75, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 3, 11, 27, 31, 73, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 7, 7, 19, 25, 105, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0,   0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 5, 5, 21, 9, 7, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 1, 15, 5, 49, 59, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 1, 1, 1, 33, 65, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 5, 15, 17, 19, 21, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 7, 11, 13, 29, 3, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 7, 5, 7, 11, 113, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,  0,   0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 5, 3, 15, 19, 61, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 3, 1, 1, 9, 27, 89, 7, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 3, 7, 31, 15, 45, 23, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
        {1, 3, 3, 9, 9, 25, 107, 39, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0,  0,   0,  0, 0, 0, 0, 0, 0, 0}};

constexpr uint64_t POLY[DIM_MAX] = {
        1,   3,   7,   11,  13,  19,  25,  37,  59,  47,  61,  55,  41,  67,
        97,  91,  109, 103, 115, 131, 193, 137, 145, 143, 241, 157, 185, 167,
        229, 171, 213, 191, 253, 203, 211, 239, 247, 285, 369, 299};

VGrid GenVGrid(uint64_t dim) {
    // Initialize V_GRID
    VGrid v_grid;

    for (uint64_t i = 0; i < DIM_MAX; i++) {
        for (uint64_t j = 0; j < LOG_MAX; j++) {
            v_grid[i][j] = V_TMPL[i][j];
        }
    }
    for (uint64_t i = 2; i < dim + 1; i++) {
        const uint64_t m = GetHighestBitPos(POLY[i - 1] >> 1);

        std::vector<uint64_t> inc(m);
        {
            uint64_t j = POLY[i - 1];
            for (uint64_t k = m; 0 < k; k--) {
                inc[k - 1] = (j != (j & (~uint64_t(1))));
                j >>= 1;
            }
        }

        for (uint64_t j = m + 1; j < LOG_MAX + 1; j++) {
            uint64_t newv = v_grid[i - 1][j - m - 1];
            uint64_t l = 1;
            for (uint64_t k = 1; k < m + 1; k++) {
                l <<= 1;
                if (inc[k - 1]) {
                    newv ^= l * v_grid[i - 1][j - k - 1];
                }
            }
            v_grid[i - 1][j - 1] = newv;
        }
    }

    uint64_t l = 1;
    for (uint64_t j = LOG_MAX - 1; 0 < j; j--) {
        l <<= 1;
        for (uint64_t k = 0; k < dim; k++) {
            v_grid[k][j - 1] *= l;
        }
    }

    return v_grid;
}

constexpr uint64_t GenRecipdDenom() {
    uint64_t l = 1;
    for (uint64_t j = LOG_MAX - 1; 0 < j; j--) {
        l <<= 1;
    }
    return 2 * l;
}

std::vector<uint64_t> GenSampleSizes2Base(
        const std::vector<uint64_t>& n_samples) {
    std::vector<uint64_t> n_samples_2base;
    n_samples_2base.reserve(n_samples.size());
    for (uint64_t i = 0; i < n_samples.size(); i++) {
        n_samples_2base.push_back(1 << GetHighestBitPos(n_samples[i]));
    }
    return n_samples_2base;
}

}  // namespace

Sobol::Sobol(uint64_t dim, uint64_t n_sample, uint64_t seed)
    : Sobol(dim, std::vector<uint64_t>(dim, n_sample), seed) {}

Sobol::Sobol(uint64_t dim, const std::vector<uint64_t>& n_samples,
             uint64_t seed)
    : DIM(dim),
      V_GRID(GenVGrid(dim)),
      RECIPD_DENOM(GenRecipdDenom()),
      N_SAMPLES(n_samples),
      N_SAMPLES_2BASE(GenSampleSizes2Base(n_samples)),
      m_seed(0),
      m_last_q(dim, 0),
      m_last_sample(dim, ~uint64_t(0)) {
    if (DIM_MAX < DIM) {
        throw std::runtime_error("Input dimension is too high");
    }
    if (N_SAMPLES_2BASE.size() != DIM) {
        throw std::runtime_error("Invalid sample size dimension");
    }
    for (uint64_t i = 0; i < N_SAMPLES_2BASE.size(); i++) {
        if (N_SAMPLES_2BASE[i] == 0) {
            throw std::runtime_error("Invalid sample size (0)");
        }
    }

    // Increment the seed
    for (uint64_t seed_tmp = 0; seed_tmp < seed; seed_tmp++) {
        nextImpl();
    }
    assert(seed == m_seed);
}

const std::vector<uint64_t>& Sobol::next() {
    // Sobol supports only 2^n. If sample is not, it should sample again
    while (true) {
        nextImpl();
        if (isValidSample()) {
            break;
        }
    }
    return m_last_sample;
}

const std::vector<uint64_t>& Sobol::getLastSample() const {
    return m_last_sample;
}

void Sobol::nextImpl() {
    const uint64_t l = GetLowestBitPos(m_seed);
    if (LOG_MAX < l) {
        // Seed became too large
        m_seed = 0;
    } else {
        m_seed++;
    }

    for (uint64_t i = 0; i < DIM; i++) {
        uint64_t s = m_last_q[i] * N_SAMPLES_2BASE[i] / RECIPD_DENOM;
        m_last_sample[i] = s;
        m_last_q[i] ^= V_GRID[i][l];
    }
}

bool Sobol::isValidSample() {
    for (uint64_t i = 0; i < DIM; i++) {
        if (N_SAMPLES[i] <= m_last_sample[i]) {
            return false;
        }
    }
    return true;
}

#endif /* end of implementation guard */
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

}  // namespace tinysobol

#endif /* end of include guard */
