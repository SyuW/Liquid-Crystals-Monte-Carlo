#ifndef SETUPHEADERDEF
#define SETUPHEADERDEF

// function declarations

// overlap functions
bool checkEllipseEllipseOverlap(double x1, double y1, double x2, double y2, double theta1, double theta2,
                                double minorAxis, double majorAxis, bool verbose);
                                
bool checkBoundaryOverlapCircle(double R, double minorAxis, double majorAxis,
                                double xc, double yc, double theta, bool debug);

// writing out
void writeOutPositionsCircle(const double majorAxis, const double minorAxis, const double boundaryRadius,
                             Matrix& posArray, const std::string fileName);

void writeOutPositionsBox(const double majorAxis, const double minorAxis,
                          const double length, const double width, Matrix& posArray, const std::string fileName);

// simulation functions
void tuneAcceptanceRate(const double rate, double& stepXY, double& stepTh);

Matrix boxHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
                                 const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
                                 Matrix posArray);

Matrix circleHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps, const double boundaryRadius,
                                    const double majorAxis, const double minorAxis, Matrix posArray);


#include "../auxiliary/list"
#include "../setup/overlap.cpp"
#include "../setup/simulators.cpp"
#include "../setup/input.cpp"
#include "../setup/output.cpp"
#include "../utils/initial.cpp"

#endif