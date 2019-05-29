/*
  License:
    MIT license.

  Notation:
    This code is based on `https://github.com/naught101/sobol_seq`.

*/

#ifndef TINY_SOBOL_H_190529
#define TINY_SOBOL_H_190529

#include <cstdint>
#include <iostream>
#include <vector>

namespace tinysobol {

// =============================================================================
// =================================== Sobol ===================================
// =============================================================================

static constexpr size_t DIM_MAX = 40;
static constexpr size_t LOG_MAX = 30;

class Sobol {
public:
    Sobol();
    Sobol(size_t dim, size_t sample_size, uint64_t seed = 0);
    Sobol(size_t dim, const std::vector<size_t>& sample_sizes,
          uint64_t seed = 0);
    ~Sobol();

    bool init(size_t dim, size_t sample_size, uint64_t seed = 0);
    bool init(size_t dim, const std::vector<size_t>& sample_sizes,
              uint64_t seed = 0);
    bool next(std::vector<uint64_t>& sample);

private:
    uint64_t m_v[DIM_MAX][LOG_MAX];
    uint64_t m_recipd_denom;
    size_t m_dim;
    std::vector<size_t> m_sample_sizes;
    std::vector<size_t> m_sample_sizes_actual;
    uint64_t m_seed;
    std::vector<uint64_t> m_lastq;
};

// *****************************************************************************
// *****************************************************************************
// **************************** Begin of Definitions ***************************
// *****************************************************************************
// *****************************************************************************
#ifdef TINY_SOBOL_IMPLEMENTATION

namespace {

// Returns the position of the high 1 bit base 2 in an integer.
uint64_t i4_bit_hi1(uint64_t n) {
    uint64_t bit = 0;
    while (n > 0) {
        bit++;
        n >>= 1;
    }
    return bit;
}

// Returns the position of the low 0 bit base 2 in an integer.
uint64_t i4_bit_lo0(uint64_t n) {
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

}  // namespace

Sobol::Sobol() {}

Sobol::Sobol(size_t dim, size_t sample_size, uint64_t seed) {
    init(dim, sample_size, seed);
}

Sobol::Sobol(size_t dim, const std::vector<size_t>& sample_sizes,
             uint64_t seed) {
    init(dim, sample_sizes, seed);
}

Sobol::~Sobol() {}

bool Sobol::init(size_t dim, size_t sample_size, uint64_t seed) {
    std::vector<size_t> sample_sizes(dim);
    for (size_t i = 0; i < dim; i++) {
        sample_sizes[i] = sample_size;
    }
    return init(dim, sample_sizes, seed);
}

bool Sobol::init(size_t dim, const std::vector<size_t>& sample_sizes,
                 uint64_t seed) {
    if (DIM_MAX < dim) {
        std::cerr << "Input dimension is too high" << std::endl;
        return false;
    }
    if (sample_sizes.size() != dim) {
        std::cerr << "Invalid sample size dimension" << std::endl;
        return false;
    }
    for (size_t i = 0; i < sample_sizes.size(); i++) {
        if (sample_sizes[i] == 0) {
            std::cerr << "Invalid sample size (0)" << std::endl;
            return false;
        }
    }

    // Initialize m_v
    for (size_t i = 0; i < DIM_MAX; i++) {
        for (size_t j = 0; j < LOG_MAX; j++) {
            m_v[i][j] = V_TMPL[i][j];
        }
    }
    for (size_t i = 2; i < dim + 1; i++) {
        const size_t m = i4_bit_hi1(POLY[i - 1] >> 1);

        std::vector<uint64_t> includ(m);
        {
            uint64_t j = POLY[i - 1];
            for (size_t k = m; 0 < k; k--) {
                includ[k - 1] = (j != (j & (~uint64_t(1))));
                j >>= 1;
            }
        }

        for (size_t j = m + 1; j < LOG_MAX + 1; j++) {
            uint64_t newv = m_v[i - 1][j - m - 1];
            uint64_t l = 1;
            for (size_t k = 1; k < m + 1; k++) {
                l <<= 1;
                if (includ[k - 1]) {
                    newv ^= l * m_v[i - 1][j - k - 1];
                }
            }
            m_v[i - 1][j - 1] = newv;
        }
    }

    uint64_t l = 1;
    for (size_t j = LOG_MAX - 1; 0 < j; j--) {
        l <<= 1;
        for (size_t k = 0; k < dim; k++) {
            m_v[k][j - 1] *= l;
        }
    }

    m_recipd_denom = 2 * l;

    m_dim = dim;
    m_sample_sizes.resize(sample_sizes.size());
    for (size_t i = 0; i < sample_sizes.size(); i++) {
        m_sample_sizes[i] = 1 << i4_bit_hi1(sample_sizes[i]);
    }
    m_sample_sizes_actual = sample_sizes;
    m_seed = seed;

    if (seed == 0) {
        m_lastq = std::vector<uint64_t>(m_dim, 0);
    } else {
        for (uint64_t seed_tmp = 0; seed_tmp < seed; seed_tmp++) {
            uint64_t l = i4_bit_lo0(seed_tmp);
            for (size_t i = 1; i < m_dim + 1; i++) {
                m_lastq[i - 1] ^= m_v[i - 1][l - 1];
            }
        }
    }

    return true;
}

bool Sobol::next(std::vector<uint64_t>& sample) {
    while (true) {
        size_t l = i4_bit_lo0(m_seed);
        if (LOG_MAX < l) {
            return false;
        }
        m_seed++;

        sample.resize(m_dim);
        bool success = true;
        for (size_t i = 1; i < m_dim + 1; i++) {
            uint64_t s =
                    m_lastq[i - 1] * m_sample_sizes[i - 1] / m_recipd_denom;
            if (m_sample_sizes_actual[i - 1] < s) {
                success = false;
            }
            sample[i - 1] = s;
            m_lastq[i - 1] ^= m_v[i - 1][l - 1];
        }
        if (success) {
            break;
        }
    }

    return true;
}

#endif /* end of implementation guard */

}  // namespace tinysobol

#endif /* end of include guard */
