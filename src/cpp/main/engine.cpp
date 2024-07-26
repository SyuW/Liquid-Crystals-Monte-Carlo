#include <iostream>
#include "setup.hpp"
#include "constants.hpp"
#include "types.hpp"

int main() 
{
    Simulation lcSim;

    double a {3};
    double b {1};

    bool test1 {checkEllipseEllipseOverlap(-b-0.1, 0, b, 0, PI/2, PI/2, b, a, true)};

    if (test1)
    {
        std::cout << "The two ellipses overlap.\n";
    }
    else
    {
        std::cout << "The two ellipses do not overlap.\n";
    }

    bool test2 {checkBoundaryOverlapCircle(4, 1, 3, 1, 0, 0, true)};

    if (test2)
    {
        std::cout << "The ellipse overlaps with boundary.\n";
    }
    else
    {
        std::cout << "The ellipse does not overlap with boundary.\n";
    }

    return 0;
}