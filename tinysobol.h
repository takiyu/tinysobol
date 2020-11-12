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
static constexpr uint32_t DIM_MAX = 40;
static constexpr uint32_t LOG_MAX = 30;
using VGrid = std::array<std::array<uint32_t, LOG_MAX>, DIM_MAX>;  // [DIM][LOG]

class Sobol {
public:
    Sobol(uint32_t dim, uint32_t n_sample, uint32_t seed = 0);
    Sobol(const std::vector<uint32_t>& n_samples, uint32_t seed = 0);

    const std::vector<uint32_t>& next();
    const std::vector<uint32_t>& getLastSample() const;

private:
    const uint32_t DIM;
    const VGrid V_GRID;
    const std::vector<uint32_t> N_SAMPLES_SHIFT;
    const std::vector<uint32_t> N_SAMPLES;

    uint32_t m_seed;
    std::vector<uint32_t> m_last_q;
    std::vector<uint32_t> m_last_sample;

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

uint32_t GetHighestBitPos(uint32_t n) {
    // Returns the position of the high 1 bit base 2 in an integer.
    uint32_t bit = 0;
    while (n > 0) {
        bit++;
        n >>= 1;
    }
    return bit;
}

uint32_t GetLowestBitPos(uint32_t n) {
    // Returns the position of the low 0 bit base 2 in an integer.
    uint32_t bit = 1;
    while (n != (n & (~uint32_t(1)))) {
        bit++;
        n >>= 1;
    }
    return bit;
}

constexpr uint32_t V_TMPL[DIM_MAX][LOG_MAX] = {
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

constexpr uint32_t POLY[DIM_MAX] = {
        1,   3,   7,   11,  13,  19,  25,  37,  59,  47,  61,  55,  41,  67,
        97,  91,  109, 103, 115, 131, 193, 137, 145, 143, 241, 157, 185, 167,
        229, 171, 213, 191, 253, 203, 211, 239, 247, 285, 369, 299};

VGrid GenVGrid(uint32_t dim) {
    // Initialize V_GRID
    VGrid v_grid;

    for (uint32_t i = 0; i < DIM_MAX; i++) {
        for (uint32_t j = 0; j < LOG_MAX; j++) {
            v_grid[i][j] = V_TMPL[i][j];
        }
    }
    for (uint32_t i = 2; i < dim + 1; i++) {
        const uint32_t m = GetHighestBitPos(POLY[i - 1] >> 1);

        std::vector<uint32_t> inc(m);
        {
            uint32_t j = POLY[i - 1];
            for (uint32_t k = m; 0 < k; k--) {
                inc[k - 1] = (j != (j & (~uint32_t(1))));
                j >>= 1;
            }
        }

        for (uint32_t j = m + 1; j < LOG_MAX + 1; j++) {
            uint32_t newv = v_grid[i - 1][j - m - 1];
            uint32_t l = 1;
            for (uint32_t k = 1; k < m + 1; k++) {
                l <<= 1;
                if (inc[k - 1]) {
                    newv ^= l * v_grid[i - 1][j - k - 1];
                }
            }
            v_grid[i - 1][j - 1] = newv;
        }
    }

    uint32_t l = 1;
    for (uint32_t j = LOG_MAX - 1; 0 < j; j--) {
        l <<= 1;
        for (uint32_t k = 0; k < dim; k++) {
            v_grid[k][j - 1] *= l;
        }
    }

    return v_grid;
}

std::vector<uint32_t> GenSampleSizesShift(
        const std::vector<uint32_t>& n_samples) {
    std::vector<uint32_t> n_samples_shift;
    n_samples_shift.reserve(n_samples.size());
    for (uint32_t i = 0; i < n_samples.size(); i++) {
        n_samples_shift.push_back(GetHighestBitPos(n_samples[i]));
    }
    return n_samples_shift;
}

}  // namespace

Sobol::Sobol(uint32_t dim, uint32_t n_sample, uint32_t seed)
    : Sobol(std::vector<uint32_t>(dim, n_sample), seed) {}

Sobol::Sobol(const std::vector<uint32_t>& n_samples,
             uint32_t seed)
    : DIM(static_cast<uint32_t>(n_samples.size())),
      V_GRID(GenVGrid(DIM)),
      N_SAMPLES(n_samples),
      N_SAMPLES_SHIFT(GenSampleSizesShift(n_samples)),
      m_seed(0),
      m_last_q(DIM, 0),
      m_last_sample(DIM, ~uint32_t(0)) {
    if (DIM_MAX < DIM) {
        throw std::runtime_error("Input dimension is too high");
    }
    for (uint32_t i = 0; i < N_SAMPLES_SHIFT.size(); i++) {
        if (N_SAMPLES_SHIFT[i] == 0) {
            throw std::runtime_error("Invalid sample size (0)");
        }
    }

    // Increment the seed
    for (uint32_t seed_tmp = 0; seed_tmp < seed; seed_tmp++) {
        nextImpl();
    }
    assert(seed == m_seed);
}

const std::vector<uint32_t>& Sobol::next() {
    // Sobol supports only 2^n. If sample is not, it should sample again
    while (true) {
        nextImpl();
        if (isValidSample()) {
            break;
        }
        std::cout << "invalid sample" << std::endl;
    }
    return m_last_sample;
}

const std::vector<uint32_t>& Sobol::getLastSample() const {
    return m_last_sample;
}

void Sobol::nextImpl() {
    const uint32_t l = GetLowestBitPos(m_seed);
    if (LOG_MAX < l) {
        // Seed became too large
        m_seed = 0;
    } else {
        m_seed++;
    }

    for (uint32_t i = 1; i < DIM + 1; i++) {
        uint32_t s = m_last_q[i - 1] >> (LOG_MAX + 1 - N_SAMPLES_SHIFT[i - 1]);
        m_last_sample[i - 1] = s;
        m_last_q[i - 1] ^= V_GRID[i - 1][l - 1];
    }
}

bool Sobol::isValidSample() {
    for (uint32_t i = 0; i < DIM; i++) {
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
