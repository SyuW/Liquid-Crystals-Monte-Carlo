#ifndef SIMULATORSHEADERDEF
#define SIMULATORSHEADERDEF

#include "../auxiliary/list"

void tuneAcceptanceRate(const double rate, double& stepXY, double& stepTh);

Matrix boxHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
                                 const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
                                 Matrix posArray);

Matrix circleHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps, const double boundaryRadius,
                                    const double majorAxis, const double minorAxis, Matrix posArray);

#endif