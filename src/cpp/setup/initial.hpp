#ifndef INITIALHEADERDEF
#define INITIALHEADERDEF

// Create initial conditions configuration for particles within box.
Matrix initializePositionsBox(const int numParticles, const double majorAxis, const double minorAxis,
                              const double height, const double width);

// CREATE initial conditions configuration for particles within circle.
Matrix initializePositionsCircle(const int numParticles, const double majorAxis, const double minorAxis,
                                 const double radius);

#endif