#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <map>
#include <vector>

#include "../utils/initial.cpp"
#include "../setup/output.cpp"
#include "simulators.cpp"


int main()
{
    // start the timer
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // system setup
    int numParticles {100};
    int numMonteCarloSteps {20000};
    double majorAxis {7};
    double minorAxis {1};

    std::string containerType { "circle" };

    std::cout << "What semi-major axis do you want?\n";
    std::cin >> majorAxis;
    std::cout << "What semi-minor axis do you want?\n";
    std::cin >> minorAxis;
    std::cout << "How many particles do you want?\n";
    std::cin >> numParticles;
    std::cout << "How many Monte Carlo steps do you want to simulate?\n";
    std::cin >> numMonteCarloSteps;

    if (containerType == "circle")
    {
        double boundaryRadius;
        std::cout << "What boundary radius do you want?\n";
        std::cin >> boundaryRadius;

        // initial positions of particles inside hard circle
        Matrix initPosArray { initializePositionsCircle(numParticles, majorAxis, minorAxis, boundaryRadius) };

        // number of particles may change if there is not enough capacity inside container
        if (numParticles != initPosArray.getNumberOfRows())
        {
            numParticles = initPosArray.getNumberOfRows();
        }

        // write out the initial positions
        writeOutPositionsCircle(majorAxis, minorAxis, boundaryRadius, initPosArray, "initialPositions.txt");
        std::cout << "Done generating/writing initial positions file.\n";

        // start the simulation
        Matrix finalPosArray = circleHardBoundaryMonteCarlo(numParticles, numMonteCarloSteps,
                                                            boundaryRadius, majorAxis, minorAxis,
                                                            initPosArray);

        // finished simulation, write out to file
        writeOutPositionsCircle(majorAxis, minorAxis, boundaryRadius, finalPosArray, "finalPositions.txt");
        std::cout << "Done writing to output file.\n";
    }

    else if (containerType == "box")
    {
        double boxHeight;
        double boxWidth;
        std::cout << "What box height do you want?\n";
        std::cin >> boxHeight;
        std::cout << "What box width do you want?\n";
        std::cin >> boxWidth;

        // initial positions of particles inside hard box
        Matrix initPosArray { initializePositionsBox(numParticles, majorAxis, minorAxis, boxHeight, boxWidth) };

        if (numParticles != initPosArray.getNumberOfRows())
        {
            numParticles = initPosArray.getNumberOfRows();
        }

        // write out the initial positions
        writeOutPositionsBox(majorAxis, minorAxis, boxHeight, boxWidth, initPosArray, "initialPositions.txt");
        std::cout << "Done generating/writing initial positions file.\n";

        // start the simulation
        Matrix finalPosArray = boxHardBoundaryMonteCarlo(numParticles, numMonteCarloSteps,
                                                         boxHeight, boxWidth, majorAxis, minorAxis,
                                                         initPosArray);

        // finished simulation, write out to file
        writeOutPositionsBox(majorAxis, minorAxis, boxHeight, boxWidth, finalPosArray, "finalPositions.txt");
    }

    // stop the clock
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}