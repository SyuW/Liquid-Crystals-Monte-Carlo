#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "src/cpp/auxiliary/list"


void writeCircleSimNotes(const double majorAxis, const double minorAxis, 
                         const double boundaryRadius,
                         const int numMonteCarloSteps,
                         const std::string fileName)
{
    std::ofstream outFile { fileName };
    outFile << "[Parameters]\n";
    outFile << "Semi-major axis: " << majorAxis << "\n";
    outFile << "Semi-minor axis: " << minorAxis << "\n";
    outFile << "Circle radius: " << boundaryRadius << "\n";
    outFile << "Number of Monte Carlo steps: " << numMonteCarloSteps << "\n";
    outFile.close();
}


void writeBoxSimNotes(const double majorAxis, const double minorAxis,
                      const double boxHeight, const double boxWidth,
                      const int numMonteCarloSteps,
                      const std::string fileName)
{
    std::ofstream outFile { fileName };
    outFile << "[Parameters]\n";
    outFile << "Semi-major axis: " << majorAxis << "\n";
    outFile << "Semi-minor axis: " << minorAxis << "\n";
    outFile << "Box height: " << boxHeight << "\n";
    outFile << "Box width: " << boxWidth << "\n";
    outFile << "Number of Monte Carlo steps: " << numMonteCarloSteps << "\n"; 
    outFile.close();
}


void writeOutPositions(Matrix& posArray, const std::string fileName)
{
    /* Write out particle positions to a txt file */

    std::ofstream outFile { fileName };
    outFile << "# x y theta\n";
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