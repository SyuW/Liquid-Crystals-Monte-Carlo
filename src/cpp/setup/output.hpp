#ifndef OUTPUTHEADERDEF
#define OUTPUTHEADERDEF

#include <string>
#include "../auxiliary/list"

void writeOutPositionsCircle(const double majorAxis, const double minorAxis, const double boundaryRadius,
                             Matrix& posArray, const std::string fileName);

void writeOutPositionsBox(const double majorAxis, const double minorAxis, const double height, const double width,
                          Matrix& posArray, const std::string fileName);

#endif