#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <map>

#include "setup.hpp"
#include "constants.hpp"
#include "types.hpp"


void writeOutPositions(Matrix& posArray, const std::string fileName)
{
    std::ofstream outFile { fileName };
    outFile << "# x y theta\n";
    // outFile << "# Semi-major axis: " << majorAxis << "\n";
    // outFile << "# Semi-minor axis: " << minorAxis << "\n";
    // outFile << "# Boundary radius: " << boundaryRadius << "\n";
    
    for (int particleIndex=0; particleIndex < posArray.getNumberOfRows(); ++particleIndex)
    {
        outFile << posArray(particleIndex, 0)
                << ","
                << posArray(particleIndex, 1)
                << ","
                << posArray(particleIndex, 2) << "\n";
    }

    outFile.close();
}


void tuneAcceptanceRate(const double rate, double& stepXY, double& stepTh)
{
    if (rate <= 0.07)
    {
        stepTh *= 0.2;
        stepXY *= 0.2;
    }
    else if (rate > 0.07 && rate <= 0.17)
    {
        stepTh *= 0.35;
        stepXY *= 0.35;
    }
    else if (rate > 0.17 && rate <= 0.27)
    {
        stepTh *= 0.6;
        stepXY *= 0.6;
    }
    else if (rate > 0.27 && rate <= 0.37)
    {
        stepTh *= 0.75;
        stepXY *= 0.75;
    }
    else if (rate > 0.37 && rate <= 0.47)
    {
        stepTh *= 0.9;
        stepXY *= 0.9;
    }
    else if (rate > 0.47 && rate <= 0.57)
    {
        stepTh *= 1;
        stepXY *= 1;
    }
    else if (rate > 0.57 && rate <= 0.67)
    {
        stepTh *= 1.1;
        stepXY *= 1.1;
    }
    else if (rate > 0.67 && rate <= 0.77)
    {
        stepTh *= 1.35;
        stepXY *= 1.35;
    }
    else if (rate > 0.77 && rate <= 0.87)
    {
        stepTh *= 1.5;
        stepXY *= 1.5;
    }
    else if (rate > 0.87 && rate <= 0.97)
    {
        stepTh *= 1.65;
        stepXY *= 1.65;
    }
    else if (rate > 0.97)
    {
        stepTh *= 1.8;
        stepXY *= 1.8;
    }
}


int main()
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // system setup
    int numParticles {100};
    int numMonteCarloSteps {20000};
    double boundaryRadius {25};
    double majorAxis {7};
    double minorAxis {1};

    std::cout << "What semi-major axis do you want?\n";
    std::cin >> majorAxis;
    std::cout << "What semi-minor axis do you want?\n";
    std::cin >> minorAxis;
    std::cout << "How many particles do you want?\n";
    std::cin >> numParticles;
    std::cout << "How many Monte Carlo steps do you want to simulate?\n";
    std::cin >> numMonteCarloSteps;

    // initial positions of particles inside container
    Matrix posArray { initializePositionsCircle(numParticles, majorAxis, minorAxis, boundaryRadius) };

    // number of particles may change if there is not enough capacity inside container
    numParticles = posArray.getNumberOfRows();

    writeOutPositions(posArray, "initialPositions.txt");

    std::cout << "Done generating/writing initial positions file.\n";

    // seed a random number generator
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> uniform_dist(-1, 1);

    double stepXY { 0.5 * boundaryRadius };
    double stepTh { PI / 2 };

    double uX {};
    double uY {};
    double uTh {};
 
    double proposedX {};
    double proposedY {};
    double proposedTh {};

    int acceptedMoves {0};
    int totalMoves {0};

    int tunePeriod { numMonteCarloSteps / 50 };

    for (int stepNo {1}; stepNo <= numMonteCarloSteps; ++stepNo)
    {
        for (int particleIndex1 {0}; particleIndex1 < numParticles; ++particleIndex1)
        {

            // tune displacements to maintain acceptance rate
            if (stepNo % tunePeriod == 0)
            {
                tuneAcceptanceRate(static_cast<double>(acceptedMoves)/totalMoves, stepXY, stepTh);
            }
            
            proposedX  = posArray(particleIndex1, 0) + stepXY * uniform_dist(rng);
            proposedY  = posArray(particleIndex1, 1) + stepXY * uniform_dist(rng);
            proposedTh = posArray(particleIndex1, 2) + stepTh * uniform_dist(rng);

            // ellipse center cannot be outside (less expensive)
            if (pow(proposedX, 2) + pow(proposedY, 2) > pow(boundaryRadius, 2))
            {
                continue;
            }

            // otherwise, check overlap with boundary (more expensive)
            else if (checkBoundaryOverlapCircle(boundaryRadius, minorAxis, majorAxis, proposedX, proposedY, proposedTh))
            {
                continue;
            }

            bool overlapVar {false};
            for (int particleIndex2 {0}; particleIndex2 < numParticles; ++particleIndex2)
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

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}