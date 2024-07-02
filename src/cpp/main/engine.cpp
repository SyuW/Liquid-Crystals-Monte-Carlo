#include <iostream>
#include "setup.hpp"
#include "constants.hpp"
#include "types.hpp"

int main() {

    double a {3};
    double b {1};

    bool test {checkOverlap(0, 0, a + b - 3, 0, 0, PI/2, b, a)};

    std::cout << test << "\n";

    return 0;
}