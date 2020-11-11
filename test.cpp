#define TINY_SOBOL_IMPLEMENTATION
#include <cassert>
#include <set>
#include <sstream>

#include "tinysobol.h"

uint64_t Pow(uint64_t v, uint64_t n) {
    if (n == 1) {
        return v;
    } else {
        return Pow(v, n - 1) * v;
    }
}

int main(int argc, char const* argv[]) {
    // Create Sobol instance
    const uint64_t DIM = 2;
    const uint64_t SAMPLE_SIZE = 4;
    tinysobol::Sobol sobol(DIM, SAMPLE_SIZE);

    std::set<std::string> histroy;

    for (uint64_t cnt = 0; cnt < Pow(SAMPLE_SIZE, DIM) + 10; cnt++) {
        // Sample
        const std::vector<uint64_t>& sample = sobol.next();

        // Value check
        std::stringstream ss;
        for (uint64_t i = 0; i < sample.size(); i++) {
            ss << sample[i] << " ";
        }
        if (histroy.count(ss.str())) {
            std::cout << "Value overlap !!!" << std::endl;
        } else {
            histroy.insert(ss.str());
        }

        // Print
        std::cout << ss.str() << std::endl;
    }

    return 0;
}
