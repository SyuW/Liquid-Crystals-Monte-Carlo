#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include "../auxiliary/list"


void writeOutPositionsCircle(const double majorAxis, const double minorAxis, const double boundaryRadius,
                             Matrix& posArray, const std::string fileName)
/*
 * Write out particle positions in circle to a txt file 
 * 
 * majorAxis
 * minorAxis
 * boundaryRadius
 * posArray
 * fileName - 
 * 
 */

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


void writeOutPositionsBox(const double majorAxis, const double minorAxis, const double height, const double width,
                          Matrix& posArray, const std::string fileName)
/*
 * Write out particle positions in box to a txt file
 *
 * majorAxis
 * minorAxis
 * height
 * width
 * posArray
 * fileName
 * 
 */
{
    std::ofstream outFile { fileName };
    outFile << "# x y theta\n";
    outFile << "# Semi-major axis: " << majorAxis << "\n";
    outFile << "# Semi-minor axis: " << minorAxis << "\n";
    outFile << "# Height: " << height << "\n";
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