#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <map>

#include "setup.hpp"
#include "constants.hpp"
#include "types.hpp"

int main() 
{
    // system setup
    int numParticlesToSimulate {100};
    double boundaryRadius {25};
    double majorAxis {5};
    double minorAxis {1};

    // initial positions of particles inside container
    Matrix posArray { initializePositionsCircle(numParticlesToSimulate, majorAxis, minorAxis, boundaryRadius) };

    std::ofstream initFile { "initPositions.txt" };
    initFile << "# x y theta\n";
    initFile << "# Semi-major axis: " << majorAxis << "\n";
    initFile << "# Semi-minor axis: " << minorAxis << "\n";
    initFile << "# Boundary radius: " << boundaryRadius << "\n";

    for (int particleIndex=0; particleIndex < posArray.getNumberOfRows(); ++particleIndex)
    {
        initFile << posArray(particleIndex, 0)
                 << ","
                 << posArray(particleIndex, 1)
                 << ","
                 << posArray(particleIndex, 2) << "\n";
    }

    std::cout << "Done generating/writing initial positions file.\n";
    initFile.close();

    // system-independent simulation parameters
    int numMonteCarloSteps {10000};

    // seed a random number generator
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> uniform_dist(-1, 1);

    double stepXY { 0.5 * boundaryRadius };
    double stepTh { PI / 2 };

    double u_x {};
    double u_y {};
    double u_t {};
 
    double proposedX { 0 };
    double proposedY { 0 };
    double proposedTh { 0 };

    int acceptedMoves {0};
    int totalMoves {0};

    for (int stepNo {1}; stepNo <= numMonteCarloSteps; ++stepNo)
    {
        for (int particleIndex1 {0}; particleIndex1 < posArray.getNumberOfRows(); ++particleIndex1)
        {
            u_x = uniform_dist(rng);
            u_y = uniform_dist(rng);
            u_t = uniform_dist(rng);

            proposedX = posArray(particleIndex1, 0) + u_x;
            proposedY = posArray(particleIndex1, 1) + u_y;
            proposedTh = posArray(particleIndex1, 2) + u_t;

            bool boundaryOverlapVar { checkBoundaryOverlapCircle(boundaryRadius, minorAxis, majorAxis, proposedX, proposedY, proposedTh) };

            if (pow(proposedX, 2) + pow(proposedY, 2) > pow(boundaryRadius, 2) || boundaryOverlapVar)
            {
                continue;
            }

            bool overlapVar {false};
            for (int particleIndex2 {0}; particleIndex2 < posArray.getNumberOfRows(); ++particleIndex2)
            {
                if (particleIndex2 == particleIndex1)
                {
                    continue;
                }

                double x2 { posArray(particleIndex2, 0) };
                double y2 { posArray(particleIndex2, 1) };
                double t2 { posArray(particleIndex2, 2) };

                if (checkEllipseEllipseOverlap(proposedX, proposedY, x2, y2, proposedTh, t2, minorAxis, majorAxis))
                {
                    overlapVar = true;
                    break;
                }
            }

            if (overlapVar)
            {
                overlapVar = false;
                continue;
            }
            else
            {
                posArray(particleIndex1, 0) = proposedX;
                posArray(particleIndex1, 1) = proposedY;
                posArray(particleIndex1, 2) = proposedTh;
            }

        }
    }

    std::cout << "Done simulation.\n";

    // now write to the outfile
    std::ofstream outFile { "finalPositions.txt" };

    if (!outFile)
    {
        std::cerr << "Could not open finalPositions.txt for writing.\n";
        return 1;
    }

    outFile << "# x y theta\n";
    outFile << "# Semi-major axis: " << majorAxis << "\n";
    outFile << "# Semi-minor axis: " << minorAxis << "\n";
    outFile << "# Boundary radius: " << boundaryRadius << "\n";
    for (int particleIndex=0; particleIndex < posArray.getNumberOfRows(); ++particleIndex)
    {
        outFile << posArray(particleIndex, 0)
                << ","
                << posArray(particleIndex, 1)
                << ","
                << posArray(particleIndex, 2) << "\n";
    }

    outFile.close();

    std::cout << "Done writing to output file.\n";

    return 0;
}