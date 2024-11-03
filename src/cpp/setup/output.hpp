#ifndef OUTPUTHEADERDEF
#define OUTPUTHEADERDEF

#include <string>
#include "src/cpp/auxiliary/list"

void writeCircleSimNotes(const double majorAxis, const double minorAxis, 
                         const double boundaryRadius,
                         const int numMonteCarloSteps,
                         const std::string fileName);

void writeBoxSimNotes(const double majorAxis, const double minorAxis,
                      const double boxHeight, const double boxWidth,
                      const int numMonteCarloSteps,
                      const std::string fileName);

void writeOutPositions(Matrix& posArray, const std::string fileName);

void writeOutPositionsCircle(const double majorAxis, const double minorAxis, const double boundaryRadius,
                             Matrix& posArray, const std::string fileName);

void writeOutPositionsBox(const double majorAxis, const double minorAxis, const double height, const double width,
                          Matrix& posArray, const std::string fileName);

#endif