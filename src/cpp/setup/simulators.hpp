#ifndef SIMULATORSHEADERDEF
#define SIMULATORSHEADERDEF

#include "src/cpp/auxiliary/list"

void tuneAcceptanceRate(const double rate, double& stepXY, double& stepTh);

Matrix boxHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
                                 const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
                                 Matrix posArray, const std::string outDir);

Matrix circleHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps, const double boundaryRadius,
                                    const double majorAxis, const double minorAxis,
                                    Matrix posArray, const std::string outDir);

#endif