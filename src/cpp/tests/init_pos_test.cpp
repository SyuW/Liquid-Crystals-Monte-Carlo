#include <chrono>
#include <iostream>
#include <fstream>

#include "src/cpp/auxiliary/list"
#include "src/cpp/setup/initial.hpp"


void test_initPositionsBox()
{
    int numParticlesToSimulate;
    double majorAxis, minorAxis;
    double boxLength, boxWidth;

    std::cout << "What length do you want for the box?\n";
    std::cin >> boxLength;

    std::cout << "What width do you want for the box?\n";
    std::cin >> boxWidth;

    std::cout << "What is the semi-major axis of the ellipse particle?\n";
    std::cin >> majorAxis;

    std::cout << "What is the semi-minor axis of the ellipse particle?\n";
    std::cin >> minorAxis;

    std::cout << "How many particles do you want to simulate?\n";
    std::cin >> numParticlesToSimulate;

    Matrix initialPositions { initializePositionsBox(numParticlesToSimulate, majorAxis, minorAxis,
                                                     boxLength, boxWidth) };
    std::ofstream outFile { "initialPositions_box.txt" };

    if (!outFile)
    {
        // print an error and exit
        std::cerr << "Could not open initialPositions_box.txt for writing.\n";
    }

    outFile << "# x y theta\n";
    for (int particleIndex=0; particleIndex < initialPositions.getNumberOfRows(); ++particleIndex)
    {
        outFile << initialPositions(particleIndex, 0)
                << ","
                << initialPositions(particleIndex, 1)
                << ","
                << initialPositions(particleIndex, 2) << "\n";
    }

    outFile.close();

    std::cout << "Initial positions file creation was successful.\n";
}


void test_initPositionsCircle()
{
    int numParticlesToSimulate {200};
    double majorAxis {3};
    double minorAxis {1};
    double boundaryRadius {25};

    std::cout << "How big do you want to make the radius of the circle boundary?\n";
    std::cin >> boundaryRadius;

    std::cout << "What is the semi-major axis of the ellipse particle?\n";
    std::cin >> majorAxis;

    std::cout << "What is the semi-minor axis of the ellipse particle?\n";
    std::cin >> minorAxis;

    std::cout << "How many particles do you want to simulate?\n";
    std::cin >> numParticlesToSimulate;

    Matrix initialPositions { initializePositionsCircle(numParticlesToSimulate, majorAxis, minorAxis, boundaryRadius) };

    std::ofstream outFile { "initialPositions_circle.txt" };

    if (!outFile)
    {
        // Print an error and exit
        std::cerr << "Could not open initialPositions_circle.txt for writing.\n";
        return;
    }

    outFile << "# x y theta\n";
    for (int particleIndex=0; particleIndex < initialPositions.getNumberOfRows(); ++particleIndex)
    {
        outFile << initialPositions(particleIndex, 0) 
                << "," 
                << initialPositions(particleIndex, 1) 
                << "," 
                << initialPositions(particleIndex, 2) << "\n";
    }

    outFile.close();

    std::cout << "Initial positions file creation was successful.\n";
}


int main()
{
    
    test_initPositionsBox();
    // test_initPositionsCircle();

    return 0;
}