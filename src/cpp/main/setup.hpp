#ifndef SETUPHEADERDEF
#define SETUPHEADERDEF

#include "src/cpp/auxiliary/list"
// #include "../setup/overlap.cpp"
#include "src/cpp/setup/simulators.hpp"
#include "src/cpp/setup/output.hpp"
#include "src/cpp/setup/initial.hpp"


// function declarations

// overlap functions
// const bool checkEllipseEllipseOverlap(const double x1, const double y1, const double x2, const double y2,
//                                       const double theta1, const double theta2,
//                                       const double minorAxis, const double majorAxis, const bool debug=false);
                                
// const bool checkBoundaryOverlapCircle(const double R, const double minorAxis, const double majorAxis,
//                                       const double xc, const double yc, const double theta, const bool debug=false);

// writing out
// void writeOutPositionsCircle(const double majorAxis, const double minorAxis, const double boundaryRadius,
//                              Matrix& posArray, const std::string fileName);

// void writeOutPositionsBox(const double majorAxis, const double minorAxis,
//                           const double length, const double width, Matrix& posArray, const std::string fileName);

// // initialization
// Matrix initializePositionsBox(const int numParticles, const double majorAxis, const double minorAxis,
//                               const double height, const double width);

// Matrix initializePositionsCircle(const int numParticles, const double majorAxis, const double minorAxis,
//                                  const double radius);

// // simulation functions
// void tuneAcceptanceRate(const double rate, double& stepXY, double& stepTh);

// Matrix boxHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
//                                  const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
//                                  Matrix posArray);

// Matrix circleHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps, const double boundaryRadius,
//                                     const double majorAxis, const double minorAxis, Matrix posArray);

#endif