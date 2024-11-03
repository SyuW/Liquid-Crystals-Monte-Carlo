#include <algorithm>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <random>
#include <string>
#include <map>
#include <vector>


#include "src/cpp/auxiliary/list"
#include "src/cpp/setup/output.hpp"
#include "src/cpp/setup/initial.hpp"
#include "src/cpp/setup/simulators.hpp"


int main()
{
    // system setup
    int numParticles {100};
    int numMonteCarloSteps {20000};
    double majorAxis {7};
    double minorAxis {1};

    // for timing
    std::chrono::steady_clock::time_point begin;

    std::vector<std::string> allowedContainerTypes;
    allowedContainerTypes.push_back("circle");
    allowedContainerTypes.push_back("box");
    std::vector<std::string> allowedBoundaryTypes;
    allowedBoundaryTypes.push_back("hard");
    allowedBoundaryTypes.push_back("periodic");

    std::string container {};
    std::cout << "What type of container do you want: ";
    bool first = true;
    for (auto const& allowed : allowedContainerTypes)
    {
        if (first) { first = false; } else { std::cout << ", "; }
        std::cout << allowed;   
    }
    std::cout << ".\n";
    std::cin >> container;

    // check that container type is allowed
    assert(std::find(allowedContainerTypes.begin(), allowedContainerTypes.end(), container)
           != allowedContainerTypes.end());

    std::cout << "What semi-major axis do you want?\n";
    std::cin >> majorAxis;
    std::cout << "What semi-minor axis do you want?\n";
    std::cin >> minorAxis;
    std::cout << "How many particles do you want?\n";
    std::cin >> numParticles;
    std::cout << "How many Monte Carlo steps do you want to simulate?\n";
    std::cin >> numMonteCarloSteps;

    if (container == "circle")
    {
        std::string outDir { "circleSimOut/" };
        std::filesystem::create_directory(outDir);
        double boundaryRadius;
        std::cout << "What boundary radius do you want?\n";
        std::cin >> boundaryRadius;
        begin = std::chrono::steady_clock::now();
        Matrix initPosArray { initializePositionsCircle(numParticles, majorAxis, minorAxis, boundaryRadius) };
        if (numParticles != initPosArray.getNumberOfRows())
        {
            numParticles = initPosArray.getNumberOfRows();
        }
        writeCircleSimNotes(majorAxis, minorAxis, boundaryRadius, numMonteCarloSteps, outDir + "circleSimNotes.txt");
        writeOutPositions(initPosArray, outDir + "initialPositions_circle.txt");
        std::cout << "Done generating/writing initial positions file.\n";
        Matrix finalPosArray = circleHardBoundaryMonteCarlo(numParticles, numMonteCarloSteps,
                                                            boundaryRadius, majorAxis, minorAxis,
                                                            initPosArray, outDir);
        writeOutPositions(finalPosArray, outDir + "finalPositions_circle.txt");
        std::cout << "Done writing final positions to output file.\n";
    }
    else if (container == "box")
    {
        // type of boundary conditions
        std::string bType {};
        std::cout << "What condition do you want to apply to the box boundary: ";
        first = true;
        for (auto const& allowed : allowedBoundaryTypes)
        {
            if (first) { first = false; } else { std::cout << ", "; }
            std::cout << allowed;
        }
        std::cout << ".\n";
        std::cin >> bType;
        // check that the container type is allowed
        assert(std::find(allowedBoundaryTypes.begin(), allowedBoundaryTypes.end(), bType)
               != allowedBoundaryTypes.end());
        // set the box height and width
        double boxHeight;
        double boxWidth;
        std::cout << "What box height do you want?\n";
        std::cin >> boxHeight;
        std::cout << "What box width do you want?\n";
        std::cin >> boxWidth;
        // create the output directory
        std::string outDir {};
        if ( bType == "hard" ) { 
            outDir = "boxHardSimOut/"; }
        else if ( bType == "periodic" ) {
            outDir = "boxPeriodicSimOut/"; }
        std::filesystem::create_directory(outDir);
        // write simulation notes file
        writeBoxSimNotes(majorAxis, minorAxis, boxHeight, boxWidth, numMonteCarloSteps, bType, outDir + "boxSimNotes.txt");
        // ready to begin the simulation, start the clock
        begin = std::chrono::steady_clock::now();
        // create the initial positions array
        Matrix initPosArray { initializePositionsBox(numParticles, majorAxis, minorAxis, boxHeight, boxWidth) };
        writeOutPositions(initPosArray, outDir + "initialPositions_box.txt");
        std::cout << "Done generating/writing initial positions file.\n";
        if (numParticles != initPosArray.getNumberOfRows()) { numParticles = initPosArray.getNumberOfRows(); }
        // start the simulation
        Matrix finalPosArray {initPosArray};
        if ( bType == "hard" ) {
            finalPosArray = boxHardBoundaryMonteCarlo(numParticles, numMonteCarloSteps,
                                                      boxHeight, boxWidth, majorAxis, minorAxis,
                                                      initPosArray, outDir);
        }
        else if ( bType == "periodic" ) {
            finalPosArray = boxPeriodicBoundaryMonteCarlo(numParticles, numMonteCarloSteps,
                                                          boxHeight, boxWidth, majorAxis, minorAxis,
                                                          initPosArray, outDir);
        }
        writeOutPositions(finalPosArray, outDir + "finalPositions_box.txt");
        std::cout << "Done writing final positions to output file.\n";
    }

    // stop the clock
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}