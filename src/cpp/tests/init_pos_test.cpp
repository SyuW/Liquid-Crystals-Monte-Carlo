#include <iostream>
#include "../utils/initial.cpp"

void testGenerate()
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

    std::ofstream outFile { "initialPositions.txt" };

    if (!outFile)
    {
        // Print an error and exit
        std::cerr << "Could not open initialPositions.txt for writing.\n";
        return;
    }

    outFile << "# x y theta\n";
    outFile << "Semi-major axis: " << majorAxis << "\n";
    outFile << "Semi-minor axis: " << minorAxis << "\n";
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

    return;
}