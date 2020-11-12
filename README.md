# TinySobol
[Sobol Sequence](https://en.wikipedia.org/wiki/Sobol_sequence) implementation with C++.

## Example

```cpp
#define TINY_SOBOL_IMPLEMENTATION
#include "tinysobol.h"

int main(int argc, char const* argv[]) {

    tinysobol::Sobol sobol({4, 4});

    for (uint32_t i = 0; i < 16; i++) {
        const std::vector<uint32_t> sample = sobol.next();
        std::cout << sample[0] << ", " << sample[1] << std::endl;
    }

    return 0;
}
```

```
$ g++ a.cpp && ./a.out
0: [0, 0]
1: [2, 2]
2: [3, 1]
3: [1, 3]
4: [1, 1]
5: [3, 3]
6: [2, 0]
7: [0, 2]
8: [0, 1]
9: [2, 3]
10: [3, 0]
11: [1, 2]
12: [1, 0]
13: [3, 2]
14: [2, 1]
15: [0, 3]
```
