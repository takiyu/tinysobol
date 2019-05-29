#include "tinysobol.h"

#include <set>
#include <sstream>
#include <cassert>

size_t pow(size_t v, size_t n) {
    if (n == 1) {
        return v;
    } else {
        return pow(v, n - 1) * v;
    }
}

int main(int argc, char const* argv[]) {

    // Create Sobol instance
    const size_t DIM = 3;
    const size_t SAMPLE_SIZE = 128;  // No overlap, if the size is pow of 2.
    tinysobol::Sobol sobol(DIM, SAMPLE_SIZE);

    std::set<std::string> histroy;

    std::vector<uint64_t> sample;
    for (size_t cnt = 0; cnt < pow(SAMPLE_SIZE, DIM); cnt++) {
        // Sample
        bool ret = sobol.next(sample);
        assert(ret);

        // Value check
        std::stringstream ss;
        for (size_t i = 0; i < sample.size(); i++) {
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
