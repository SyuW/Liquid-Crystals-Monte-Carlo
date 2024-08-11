#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include "../auxiliary/list"


void writeOutPositionsCircle(const double majorAxis, const double minorAxis, const double boundaryRadius,
                             Matrix& posArray, const std::string fileName)
{
    std::ofstream outFile { fileName };
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
}


void writeOutPositionsBox(const double majorAxis, const double minorAxis, const double length, const double width,
                          Matrix& posArray, const std::string fileName)
{
    std::ofstream outFile { fileName };
    outFile << "# x y theta\n";
    outFile << "# Semi-major axis: " << majorAxis << "\n";
    outFile << "# Semi-minor axis: " << minorAxis << "\n";
    outFile << "# Length: " << length << "\n";
    outFile << "# Width: " << width << "\n";

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