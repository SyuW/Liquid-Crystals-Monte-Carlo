#ifndef INITIALHEADERDEF
#define INITIALHEADERDEF

Matrix initializePositionsBox(const int numParticles, const double majorAxis, const double minorAxis,
                              const double height, const double width);

Matrix initializePositionsCircle(const int numParticles, const double majorAxis, const double minorAxis,
                                 const double radius);

#endif