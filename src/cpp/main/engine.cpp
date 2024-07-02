#include <iostream>
#include "setup.hpp"
#include "constants.hpp"
#include "types.hpp"

int main() {

    double a {3};
    double b {1};

    bool test {checkOverlap(-b, 0, b+0.01, 0, PI/2, PI/2, b, a, true)};

    if (test)
    {
        std::cout << "The two ellipses overlap.\n";
    }
    else
    {
        std::cout << "The two ellipses do not overlap.\n";
    }

    std::cout << test << "\n";

    return 0;
}