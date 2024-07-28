#include <iomanip>
#include <iostream>
#include <random>
#include <map>

#include "setup.hpp"
#include "constants.hpp"
#include "types.hpp"

int main() 
{
    int numParticlesToSimulate {3};
    double boundaryRadius {25};
    double majorAxis {3};
    double minorAxis {1};

    Matrix initPos { initializePositionsCircle(numParticlesToSimulate, majorAxis, minorAxis, boundaryRadius) };

    // seed a random number generator
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> uniform_dist(0, 1);

    std::map<int, int> hist;
    for (int n = 0; n != 100000; ++n)
        ++hist[std::round(uniform_dist(rng))];
 
    std::cout << "Uniform distribution between 0 and 10.\n"
              << std::fixed << std::setprecision(1);
    for (auto [x, y] : hist)
        std::cout << std::setw(2) << x << ' ' << std::string(y / 200, '*') << '\n';

    return 0;
}